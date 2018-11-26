/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"
#include <numeric>

namespace SSAGES
{
	//! Collective variable on the distance between two particles' centers of mass.
	/*!
	 *  Collective variable on two particle positions. This CV will return the
	 *  distance between two specific atom groups of the simulation.
	 *
	 *  \ingroup CVs
	 */
	class ParticleSeparationCV : public CollectiveVariable
	{
	private:
		Label group1_; //!< Atom ID's of group 1.
		Label group2_; //!< Atom ID's of group 2.
		
		//! Each dimension determines if it is applied by the CV.
		Bool3 dim_; 

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
		group1_(group1), group2_(group2), dim_{true, true, true}
		{}

		//! Constructor
		/*!
		 * \param group1 Vector of atom ID's of the first particle.
		 * \param group2 Vector of atom ID's of the second particle.
		 * \param dimx If \c True, include x dimension.
		 * \param dimy If \c True, include y dimension.
		 * \param dimz If \c True, include z dimension.
		 *
		 * Construct a particle separation CV.
		 *
		 * \todo Bounds needs to be an input.
		 */
		ParticleSeparationCV(const Label& group1, const Label& group2, bool dimx, bool dimy, bool dimz) :
		group1_(group1), group2_(group2), dim_{dimx, dimy, dimz}
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
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});
			boxgrad_ = Matrix3::Zero();

			// Get centers of mass.
			auto mtot1 = snapshot.TotalMass(idx1);		
			auto mtot2 = snapshot.TotalMass(idx2);	
			Vector3 com1 = snapshot.CenterOfMass(idx1, mtot1);
			Vector3 com2 = snapshot.CenterOfMass(idx2, mtot2);

			// Account for pbc. 
			Vector3 rij = snapshot.ApplyMinimumImage(com1 - com2).cwiseProduct(dim_.cast<double>());

			// Compute gradient.
			val_ = rij.norm();

			for(auto& id : idx1)
			{
				grad_[id] = rij/val_*masses[id]/mtot1;
				boxgrad_ += grad_[id]*rij.transpose();
			}
			
			for(auto& id : idx2)
				grad_[id] = -rij/val_*masses[id]/mtot2;
		}
		
		//! \copydoc CollectiveVariable::BuildCV()
		static ParticleSeparationCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::ParticleSeparationCV.c_str(),
			              JsonSchema::ParticleSeparationCV.c_str() + JsonSchema::ParticleSeparationCV.size(),
			              &schema, NULL);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<int> group1, group2;
			
			for(auto& s : json["group1"])
				group1.push_back(s.asInt());

			for(auto& s : json["group2"])
				group2.push_back(s.asInt());

			ParticleSeparationCV* c;
			if(json.isMember("dimension"))
			{
				auto fixx = json["dimension"][0].asBool();
				auto fixy = json["dimension"][1].asBool();
				auto fixz = json["dimension"][2].asBool();
				c = new ParticleSeparationCV(group1, group2, fixx, fixy, fixz);
			}
			else
				c = new ParticleSeparationCV(group1, group2);
			
			return c;
		}
	};
}
