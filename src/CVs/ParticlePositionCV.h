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
		Label atomids_; 

		//! Target point in space.
		Vector3 point_;

		//! Each dimension determines if a constraint is applied by the CV.
		Bool3 dim_;
	
	public:
		//! Constructor
		/*!
		 * \param atomids Vector of atom ID's.
		 * \param position Point in (1,2,3)D space for the distance calculation.
		 *
		 * Construct a particle position CV.
		 *
		 * \todo Bounds needs to be an input.
		 */
		ParticlePositionCV(const Label& atomids, const Vector3& position) :
		atomids_(atomids), point_(position), dim_{true, true, true}
		{
		}

		//! Constructor
		/*!
		 * \param atomids Vector of atom ID's.
		 * \param position Point in (1,2,3)D space for the distance calculation.
		 * \param dimx If \c True, include x dimension.
		 * \param dimy If \c True, include y dimension.
		 * \param dimz If \c True, include z dimension.
		 *
		 * Construct a particle position CV.
		 *
		 * \todo Bounds needs to be an input.
		 */
		ParticlePositionCV(const Label& atomids, const Vector3& position, bool dimx, bool dimy, bool dimz) :
		atomids_(atomids), point_(position), dim_{dimx, dimy, dimz}
		{
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			auto n = atomids_.size();

			// Make sure atom ID's are on at least one processor. 
			std::vector<int> found(n, 0);
			for(size_t i = 0; i < n; ++i)
			{
				if(snapshot.GetLocalIndex(atomids_[i]) != -1)
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
			snapshot.GetLocalIndices(atomids_, &idx);

			// Get data from snapshot. 
			auto n = snapshot.GetNumAtoms();
			const auto& masses = snapshot.GetMasses();

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});

			// Compute total and center of mass.
			auto masstot = snapshot.TotalMass(idx);
			Vector3 com = snapshot.CenterOfMass(idx);

			// Compute difference between point and account for requested 
			// dimensions only.
			Vector3 dx = snapshot.ApplyMinimumImage(com - point_).cwiseProduct(dim_.cast<double>());
			val_ = dx.norm();

			// If distance is zero, we have nothing to do.
			if(val_ == 0)
				return;

			for(auto& id : idx)
				grad_[id] = dx/val_*masses[id]/masstot;
		}
		
		//! \copydoc CollectiveVariable::BuildCV()
		static ParticlePositionCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::ParticlePositionCV.c_str(),
			              JsonSchema::ParticlePositionCV.c_str() + JsonSchema::ParticlePositionCV.size(),
			              &schema, nullptr);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			Label atomids;
			for(auto& id : json["atom_ids"])
				atomids.push_back(id.asInt());
			
			Vector3 position;
			position[0] = json["position"][0].asDouble();
			position[1] = json["position"][1].asDouble();
			position[2] = json["position"][2].asDouble();

			ParticlePositionCV* c;

			if(json.isMember("dimension"))
			{
				auto dimx = json["dimension"][0].asBool();
				auto dimy = json["dimension"][1].asBool();
				auto dimz = json["dimension"][2].asBool();
				c = new ParticlePositionCV(atomids, position, dimx, dimy, dimz);
			}
			else
				c = new ParticlePositionCV(atomids, position);

			return c;
		}
	};
}
