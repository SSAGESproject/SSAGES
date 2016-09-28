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

#include "CollectiveVariable.h"
#include "Drivers/DriverException.h" 


namespace SSAGES
{
	//! Collective variable on the distance between two particles' centers of mass.
	/*!
	 *  Collective variable on two particle positions. This CV will return the
	 *  distance between two specific atom groups of the simulation.
	 *
	 *  \ingroup CVs
	 */
	class  ParticleSeparationCV : public CollectiveVariable
	{
	private:
		Label group1_; //! Atom ID's of group 1. 
		Label group2_; //! Atom ID's of group 2.

	public:
		//! Constructor
		/*!
		 * \param group1 Vector of atom ID's of the first particle.
		 * \param group2 Vector of atom ID's of the second particle.
		 *
		 * Construct a paarticle separation CV.
		 *
		 * \todo bounds needs to be an input.
		 */
		ParticleSeparationCV(const Label& group1, const Label& group2) : 
		group1_(group1), group2_(group2)
		{}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			auto n1 = group1_.size(), n2 = group2_.size();

			// Make sure atom ID's are on at least one processor. 
			std::vector<int> found1(n1, 0), found2(n2, 0);
			for(size_t i = 0; i < n1; ++i)
			{
				if(snapshot.GetLocalIndex(group1_[i]) != -1)
					found1[i] = 1;
			}
			
			for(size_t i = 0; i < n2; ++i)
			{
				if(snapshot.GetLocalIndex(group2_[i]) != -1)
					found2[i] = 1;
			}

			MPI_Allreduce(MPI_IN_PLACE, found1.data(), n1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
			MPI_Allreduce(MPI_IN_PLACE, found2.data(), n2, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
	
			unsigned ntot1 = std::accumulate(found1.begin(), found1.end(), 0, std::plus<int>());
			unsigned ntot2 = std::accumulate(found2.begin(), found2.end(), 0, std::plus<int>());
			if(ntot1 != n1)
				throw BuildException({
					"ParticleSeparationCV: Expected to find " + 
					to_string(n1) + 
					" atoms in particle 1, but only found " + 
					to_string(ntot1) + "."
				});			

			if(ntot2 != n2)
				throw BuildException({
					"ParticleSeparationCV: Expected to find " + 
					to_string(n2) + 
					" atoms in particle 2, but only found " + 
					to_string(ntot2) + "."
				});			
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Get local atom indices.
			std::vector<int> idx1, idx2;
			snapshot.GetLocalIndices(group1_, &idx1);
			snapshot.GetLocalIndices(group2_, &idx2);

			// Get data from snapshot. 
			auto n = snapshot.GetNumAtoms();
			const auto& masses = snapshot.GetMasses();

			// Initialize gradient.
			std::fill(_grad.begin(), _grad.end(), Vector3{0,0,0});
			_grad.resize(n, Vector3{0,0,0});

			// Get centers of mass.
			auto mtot1 = snapshot.TotalMass(idx1);		
			auto mtot2 = snapshot.TotalMass(idx2);	
			Vector3 com1 = snapshot.CenterOfMass(idx1, mtot1);
			Vector3 com2 = snapshot.CenterOfMass(idx2, mtot2);

			// Account for pbc. 
			Vector3 rij = snapshot.ApplyMinimumImage(com1 - com2);

			// Compute gradient.
			_val = rij.norm();

			for(auto& id : idx1)
				_grad[id] = rij/_val*masses[id]/mtot1;
			
			for(auto& id : idx2)
				_grad[id] = -rij/_val*masses[id]/mtot2;
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "ParticleSeparation";
			
			for(auto& id : group1_)
				json["group1"].append(id);

			for(auto& id : group2_)
				json["group2"].append(id);

			for(auto& bound : _bounds)
				json["bounds"].append(bound);
		}
	};
}