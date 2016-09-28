/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
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

#include "Drivers/DriverException.h"
#include "CollectiveVariable.h"

namespace SSAGES
{
	//! Collective variable on an particle position.
	/*!
	 * This CV will return the distance of a particle (collection of atoms) 
	 * from a particular point in (1,2,3)-dimensional space.
	 *
	 * \ingroup CVs
	 */
	class ParticlePositionCV : public CollectiveVariable
	{
	private:
		//! Vector of atom ids of interest.
		Label _atomids; 

		//! Target point in space.
		Vector3 _point;

		//! Each dimension determines if a constraint is applied by the CV.
		Bool3 _fix; 
	
	public:
		//! Constructor
		/*!
		 * \param atomids Vector of atom ID's.
		 * \param position Point in (1,2,3)D space for the distance calculation.
		 * \param fixx If \c False, constrain x dimension.
		 * \param fixy If \c False, constrain y dimension.
		 * \param fixz If \c False, constrain z dimension.
		 *
		 * Construct a particle position CV. If a dimension is constrained, this
		 * dimension does not contribute to the distance calculation. For
		 * example, if the y- and z-dimension are constrained, the CV calculates
		 * only the distance in x-direction.
		 *
		 * \todo Bounds needs to be an input.
		 */
		ParticlePositionCV(const Label& atomids, const Vector3& position, bool fixx, bool fixy, bool fixz) : 
		_atomids(atomids), _point(position), _fix{fixx, fixy, fixz}
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
					"ParticlePositionCV: Expected to find " + 
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

			// Initialize gradient.
			std::fill(_grad.begin(), _grad.end(), Vector3{0,0,0});
			_grad.resize(n, Vector3{0,0,0});

			// Compute total and center of mass.
			auto masstot = snapshot.TotalMass(idx);
			Vector3 com = snapshot.CenterOfMass(idx, masstot);

			// Compute difference between point and account for requested 
			// dimensions only.
			Vector3 dx = snapshot.ApplyMinimumImage(com - _point).cwiseProduct(_fix.cast<double>());
			_val = dx.norm();

			// If distance is zero, we have nothing to do.
			if(_val == 0)
				return;

			for(auto& id : idx)
				_grad[id] = dx/_val*masses[id]/masstot;
		}
	
		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "ParticlePosition";
			json["fix"][0] = _fix[0];
			json["fix"][1] = _fix[1];
			json["fix"][2] = _fix[2];
			
			for(auto& id : _atomids)
				json["atom_ids"].append(id);

			for(auto& bound : _bounds)
				json["bounds"].append(bound);			
		}
	};
}