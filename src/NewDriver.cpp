/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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
#include "CVs/CVManager.h"
#include "Drivers/DriverException.h"
#include "Validator/ObjectRequirement.h"
#include "Methods/Method.h"
#include "NewDriver.h"
#include "schema.h"
#include <mxx/comm.hpp>

using namespace Json;

namespace SSAGES
{
	NewDriver::NewDriver(
		const MPI_Comm& world, const MPI_Comm& comm, uint walkerid,
		const std::vector<Method*>& methods, CVManager* cvmanager) : 
	world_(world), comm_(comm), walkerid_(walkerid), methods_(methods), cvmanager_(cvmanager)
	{	
	}

	NewDriver* NewDriver::Build(const Value& json, const MPI_Comm& world)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;

		// Parse and validate top level schema. This just 
		// makes sure the proper fields exist and the correct 
		// types are specified in the input files.
		reader.parse(JsonSchema::Simulation, schema);
		validator.Parse(schema, "#");
		validator.Validate(json, "#");
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());
		
		// Get number of desired walkers and create array of input files.
		auto nwalkers = json.get("walkers", 1).asInt();
		if(json["input"].isArray() && json["input"].size() != nwalkers)
			throw BuildException({"#/input: Number of inputs do not match requested walkers."});
		
		std::vector<std::string> inputs;
		for(int i = 0; i < nwalkers; ++i)
		{
			if(json["input"].isArray())
				inputs.push_back(json["input"][i].asString());
			else
				inputs.push_back(json["input"].asString());
		}

		// Get basic MPI info and verify that the total number of 
		// cores are divisble by number of walkers. 
		auto wcomm = mxx::comm(world);
		if(wcomm.size() % nwalkers != 0) 
			throw BuildException({"#/walkers: Allocated processors not divisible by number of walkers."});
		
		// Split communicators.
		int walkerid = wcomm.rank()/nwalkers;
		auto comm = wcomm.split(walkerid);
	}
}