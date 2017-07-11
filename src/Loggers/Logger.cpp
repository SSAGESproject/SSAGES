/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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

#include <stdexcept>
#include "Logger.h"
#include "Snapshot.h"
#include "CVs/CVManager.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "schema.h"

using namespace Json;

namespace SSAGES
{
	void Logger::PreSimulation(Snapshot* snapshot, const CVManager& cvmanager)
	{
		if(comm_.rank() == 0)
		{
			if(append_)
		 		log_.open(filename_.c_str(), std::ofstream::out | std::ofstream::app);	
			else
			{
				// Write out header.
				log_.open(filename_.c_str(), std::ofstream::out);
				log_ << "#"; 
				log_ << "Iteration "; 

				auto cvs = cvmanager.GetCVs(cvmask_); 
				for(size_t i = 0; i < cvs.size() - 1; ++i)
					log_ << "cv_" + std::to_string(i) << " ";
				log_ << "cv_" + std::to_string(cvs.size() - 1) << std::endl;
			}
		}
	}

	void Logger::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		// Get necessary info. 
		auto cvs = cvmanager.GetCVs(cvmask_);
		if(comm_.rank() ==0)
		{
			log_.precision(8);
			log_ << snapshot->GetIteration() << " ";

			// Print out CV values.
			for(size_t i = 0; i < cvs.size() - 1; ++i)
				log_ << cvs[i]->GetValue() << " ";
			log_ << cvs.back()->GetValue() << std::endl;
		}
	}

	void Logger::PostSimulation(Snapshot* snapshot, const class CVManager& cvmanager)
	{
	}

	Logger* Logger::Build(const Value& json,
	                      const MPI_Comm& world,
	                      const MPI_Comm& comm,
	                      const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;

		reader.parse(JsonSchema::Logger, schema);
		validator.Parse(schema, path);
		
		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());
		
		auto freq = json.get("frequency", 1).asInt();
		std::string name = "cvlog.dat";
		
		//TODO walker id should be obtainable in method as
		//     opposed to calculated like this. 
		bool ismulti = mxx::comm(world).size() > mxx::comm(comm).size(); 
		uint wid = mxx::comm(world).rank()/mxx::comm(comm).size();
		if(json["output_file"].isArray())
			name = json["output_file"][wid].asString(); 
		else if(ismulti)
			throw std::invalid_argument(path + ": Multi-walker simulations require a separate output file for each.");
		else
			name = json["output_file"].asString();

		auto* l = new Logger(freq, name, world, comm);
		l->SetAppend(json.get("append", false).asBool());

		// Load cv mask. 
		std::vector<uint> cvmask; 
		for(auto& v : json["cvs"])
		{
			if(v.isString())
			{
				auto id = CVManager::LookupCV(v.asString());
				if(id == -1)
					throw std::invalid_argument(path + ": CV mask name \"" + v.asString() + "\" does not exist.");
				
				cvmask.push_back(CVManager::LookupCV(v.asString()));
			}
			else if(v.isIntegral() && v.asInt() >= 0)
				cvmask.push_back(v.asUInt());
			else
				throw std::invalid_argument(path + ": CV mask must contain strings or unsigned integers.");
		}
		l->SetCVMask(cvmask);

		return l;
	}

}
