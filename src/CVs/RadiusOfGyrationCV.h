/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Yamil Colon <yamilcolon2015@u.northwestern.edu>
 *                Hythem Sidky <hsidky@nd.edu>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once 

#include "CollectiveVariable.h"
#include "Drivers/DriverException.h"
#include <cmath>

namespace SSAGES
{
	//! Collective variable on the radius of gyration.
	/*!
	 * Collective variable on radius of gyration. This will return the
	 * Rg for a group of atoms.
	 *
	 * \ingroup CVs
	 */
	class RadiusOfGyrationCV : public CollectiveVariable
	{
	private:
		Label _atomids; //!< IDs of the atoms used for Rg calculation

	public:
		//! Constructor.
		/*!
		 * \param atomids IDs of the atoms defining Rg.
		 *
		 * Construct a radius of gyration CV.
		 *
		 * \todo Bounds needs to be an input
		 */
		RadiusOfGyrationCV(const Label& atomids) :
		_atomids(atomids)
		{
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			auto n = _atomids.size();

			// Make sure atom ID's are on at least one processor. 
			std::vector<int> found(n, 0);
			for(size_t i = 0; i < n; ++i)
			{
				if(snapshot.GetLocalIndex(_atomids[i]) != -1)
					found[i] = 1;
			}

			MPI_Allreduce(MPI_IN_PLACE, found.data(), n, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
			unsigned ntot = std::accumulate(found.begin(), found.end(), 0, std::plus<int>());
			if(ntot != n)
				throw BuildException({
					"RadiusOfGyrationCV: Expected to find " + 
					to_string(n) + 
					" atoms, but only found " + 
					to_string(ntot) + "."
				});			
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Get local atom indices and compute COM. 
			std::vector<int> idx;
			snapshot.GetLocalIndices(_atomids, &idx);

			// Get data from snapshot. 
			auto n = snapshot.GetNumAtoms();
			const auto& masses = snapshot.GetMasses();
			const auto& pos = snapshot.GetPositions();

			// Initialize gradient.
			_grad.resize(n, Vector3{0,0,0});
			
			// Compute total and center of mass.
			auto masstot = snapshot.TotalMass(idx);
			Vector3 com = snapshot.CenterOfMass(idx, masstot);		
			auto mrsq = 0.;
			for(auto& i : idx)
			{
				Vector3 ricm = snapshot.ApplyMinimumImage(pos[i] - com);
				mrsq += masses[i]*ricm.squaredNorm();

				// Pre-compute part of gradient for efficiency.
				auto mfrac = masses[i]/masstot;
				_grad[i] = mfrac*ricm*(1.-mfrac);
			}

			MPI_Allreduce(MPI_IN_PLACE, &mrsq, 1, MPI_DOUBLE, MPI_SUM, snapshot.GetCommunicator());
			_val = std::sqrt(mrsq/masstot);

			for(auto& i : idx)
				_grad[i] /= _val;
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{
			json["type"] = "RadiusOfGyration";
			
			for(auto& id : _atomids)
				json["atom_ids"].append(id);

			for(auto& bound : _bounds)
				json["bounds"].append(bound);		}
	};
}
