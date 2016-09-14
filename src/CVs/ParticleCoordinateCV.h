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
	//! Collective variable on a particle coordinate. 
	/*!
	 * This will return the value of either the x, y, or z coordinate, depending
	 * on the user specification for a defined particle, which is a collection of 
	 * one or more atoms.
	 *
	 * \ingroup CVs
	 */
	class ParticleCoordinateCV : public CollectiveVariable
	{
	private:
		//! IDs of atoms of interest. 
	 	Label _atomids; 

	 	//! Index of dimnesion.
	 	Dimension _dim;

	public:
		//! Constructor.
		/*!
		 * \param atomids Atom ID's of interest.
		 * \param index Index of dimension.
		 *
		 * Construct a particle coordinate CV. The atomids specify a vector of the atom
		 * ID's of interest, and index specifies the dimension to report (x,y,z).
		 *
		 * \todo Bounds needs to be an input.
		 */	 	
		ParticleCoordinateCV(const Label& atomids, Dimension dim) : 
		_atomids(atomids), _dim(dim)
		{}

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
					"ParticleCoordinateCV: Expected to find " + 
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

			// Assign CV value. 
			switch(_dim)
			{
				case Dimension::x:
					_val = com[0];
					break;
				case Dimension::y:
					_val = com[1];
					break;
				case Dimension::z:
					_val = com[2];
					break;
			}

			// Assign gradient to appropriate atoms.
			for(auto& id : idx)
			{
				_grad[id][0] = (_dim == Dimension::x) ? masses[id]/masstot : 0;
				_grad[id][1] = (_dim == Dimension::y) ? masses[id]/masstot : 0;
				_grad[id][2] = (_dim == Dimension::z) ? masses[id]/masstot : 0;
			}
		}

		//! Return value taking periodic boundary conditions into account.
		/*!
		 * \param Location Get wrapped value of this location.
		 *
		 * \return Input value
		 *
		 * The AtomCoordinate CV does not consider periodic boundary
		 * conditions. Thus, this function always returns the input value.
		 */
		double GetPeriodicValue(double Location) const override
		{
			return Location;
		}

		//! Return difference considering periodic boundaries.
		/*!
		 * \param Location Calculate difference of CV value to this location.
		 *
		 * \return Direct difference.
		 *
		 * As the AtomCoordinate CV does not consider periodic boundary
		 * conditions, the difference between the current value of the CV and
		 * another value is always the direct difference.
		 */
		double GetDifference(const double Location) const override
		{
			return _val - Location;
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "ParticleCoordinate";			
			switch(_dim)
			{
				case Dimension::x:
					json["dimension"] = "x";
					break;
				case Dimension::y:
					json["dimension"] = "y";
					break;
				case Dimension::z:
					json["dimension"] = "z";
					break;
			}
			
			for(auto& id : _atomids)
				json["atom_ids"].append(id);

			for(auto& bound : _bounds)
				json["bounds"].append(bound);
		}

	 };
}