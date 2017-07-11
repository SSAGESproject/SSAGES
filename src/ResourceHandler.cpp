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
#include <mxx/comm.hpp>
#include "CVs/CVManager.h"
#include "Drivers/DriverException.h"
#include "Loggers/Logger.h"
#include "Validator/ObjectRequirement.h"
#include "Methods/Method.h"
#include "ResourceHandler.h"
#include "Snapshot.h"
#include "Hook.h"
#include "schema.h"

using namespace Json;

namespace SSAGES
{
	ResourceHandler::ResourceHandler(
		mxx::comm&& world, mxx::comm&& comm, uint walkerid,
		const std::vector<Method*>& methods, CVManager* cvmanager) : 
	world_(std::move(world)), comm_(std::move(comm)), walkerid_(walkerid), nwalkers_(1), 
	methods_(methods), cvmanager_(cvmanager), hook_(nullptr), inputs_(0)
	{
		snapshot_ = new Snapshot(comm_, walkerid);
	}

	void ResourceHandler::ConfigureHook(class Hook* hook)
	{
		hook->SetSnapshot(snapshot_);
		hook->SetCVManager(cvmanager_);
		for(auto& m : methods_)
			hook->AddListener(m);
		if(logger_ != nullptr) hook->AddListener(logger_);
		hook_ = hook;
	}

	ResourceHandler* ResourceHandler::Build(const Value& json, const MPI_Comm& world)
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
		int ppw = wcomm.size()/nwalkers;
		int walkerid = wcomm.rank()/ppw;
		auto comm = wcomm.split(walkerid);

		// Build collective variables. 
		auto* cvmanager = new CVManager();
		int icv = 0;
		for(auto& cv : json["CVs"])
		{
			// Register name with CV manager if it exists.
			if(cv.isMember("name"))
				CVManager::AddCVtoMap(cv["name"].asString(), icv);

			// Build CV and add to manager.
			cvmanager->AddCV(CollectiveVariable::BuildCV(cv, "#/CVs"));
			++icv;
		}

		// Build methods. 
		std::vector<Method*> methods; 
		for(auto& m : json["methods"])
			methods.push_back(Method::BuildMethod(m, world, comm, "#/methods"));
		
		// Build logger. 
		Logger* l = nullptr; 
		if(json.isMember("logger"))
			l = Logger::Build(json["logger"], world, comm, "#/logger");

		auto* rh = new ResourceHandler(std::move(world), std::move(comm), walkerid, methods, cvmanager);
		rh->inputs_ = inputs;
		rh->nwalkers_ = nwalkers;
		rh->logger_ = l;
		return rh;
	}

	ResourceHandler::~ResourceHandler()
	{
		delete snapshot_; 
		delete cvmanager_;
		for(auto& m : methods_)
			delete m; 
	}
}