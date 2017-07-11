/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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
#include "Umbrella.h"
#include "Snapshot.h" 
#include "CVs/CVManager.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "schema.h"
#include <iostream>

using namespace Json;

namespace SSAGES
{
	void Umbrella::PreSimulation(Snapshot* /* snapshot */, const CVManager& cvmanager)
	{
		if(comm_.rank() == 0)
		{
			if(append_)
		 		umbrella_.open(filename_.c_str(), std::ofstream::out | std::ofstream::app);
			else
			{
				// Write out header.
				umbrella_.open(filename_.c_str(), std::ofstream::out);
				umbrella_ << "#"; 
				umbrella_ << "Iteration "; 

				auto cvs = cvmanager.GetCVs(cvmask_); 
				for(size_t i = 0; i < cvs.size(); ++i)
					umbrella_ << "cv_" + std::to_string(i) << " ";

				for(size_t i = 0; i < cvs.size() - 1; ++i)
					umbrella_ << "center_" + std::to_string(i) << " ";  
				umbrella_ << "center_" + std::to_string(cvs.size() - 1) << std::endl;
			}
		}
	}

	void Umbrella::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		// Get necessary info. 
		auto cvs = cvmanager.GetCVs(cvmask_);
		auto& forces = snapshot->GetForces();
		auto& virial = snapshot->GetVirial();
		auto& H = snapshot->GetHMatrix();

		for(size_t i = 0; i < cvs.size(); ++i)
		{
			// Get current CV and gradient.
			auto& cv = cvs[i];
			auto& grad = cv->GetGradient();
			auto& boxgrad = cv->GetBoxGradient();
			// Compute dV/dCV.
			auto center = GetCurrentCenter(snapshot->GetIteration(), i);
			auto D = kspring_[i]*cv->GetDifference(center);

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
				forces[j] -= D*grad[j];
			
			// Update virial. 
			virial += D*boxgrad;
		}

		if(snapshot->GetIteration() % outfreq_ == 0)
			PrintUmbrella(cvs, snapshot->GetIteration());
	}

	void Umbrella::PostSimulation(Snapshot*, const CVManager&)
	{
		if(comm_.rank() ==0)
			umbrella_.close();
	}

	void Umbrella::PrintUmbrella(const CVList& cvs, uint iteration)
	{
		if(comm_.rank() ==0)
		{
			umbrella_.precision(8);
			umbrella_ << iteration << " ";

			// Print out CV values first.
			for(auto& cv : cvs)
				umbrella_ << cv->GetValue() << " ";
			
			// Print out target (center) of each CV. 
			for(size_t i = 0; i < cvs.size() - 1; ++i)
				umbrella_ << GetCurrentCenter(iteration, i) << " "; 
			umbrella_ << GetCurrentCenter(iteration, cvs.size() - 1);

			umbrella_ << std::endl;
		}
	}

	Umbrella* Umbrella::Build(const Json::Value& json, 
			    		          const MPI_Comm& world,
					              const MPI_Comm& comm,
					              const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;

		reader.parse(JsonSchema::UmbrellaMethod, schema);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		std::vector<double> ksprings;
		for(auto& s : json["ksprings"])
			ksprings.push_back(s.asDouble());

		std::vector<double> centers0, centers1;
		auto timesteps = json.get("timesteps", 0).asInt();
		if(json.isMember("centers"))
		{
			for(auto& s : json["centers"])
				centers0.push_back(s.asDouble());
		}
		else if(json.isMember("centers0") && json.isMember("centers1") && json.isMember("timesteps"))
		{
			for(auto& s : json["centers0"])
				centers0.push_back(s.asDouble());

			for(auto& s : json["centers1"])
				centers1.push_back(s.asDouble());
		}
		else
			throw BuildException({"Either \"centers\" or \"timesteps\", \"centers0\" and \"centers1\" must be defined for umbrella."});

		if(ksprings.size() != centers0.size())
			throw BuildException({"Need to define a spring for every center or a center for every spring!"});

		auto freq = json.get("frequency", 1).asInt();
		std::string name = "umbrella.dat";

		//TODO walker id should be obtainable in method as
		//     opposed to calculated like this. 
		uint wid = mxx::comm(world).rank()/mxx::comm(comm).size(); 
		if(json["output_file"].isArray())
			name = json["output_file"][wid].asString(); 
		else
			name = json["output_file"].asString();

		Umbrella* m = nullptr;
		if(timesteps == 0)
			m = new Umbrella(world, comm, ksprings, centers0, name, freq);
		else
			m = new Umbrella(world, comm, ksprings, centers0, centers1, timesteps, name, freq);

		m->SetOutputFrequency(json.get("output_frequency",0).asInt());
		m->SetAppend(json.get("append", false).asBool());
		
		return m;
	}
}
