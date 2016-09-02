/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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
#include "SimObservable.h"
#include "SimObserver.h"
#include "json/json.h"
#include "schema.h"
#include "../Validator/ObjectRequirement.h"
#include "../Validator/ArrayRequirement.h"
#include "../Observers/JSONObserver.h"
#include "../Utility/BuildException.h"

using namespace Json;

namespace SSAGES
{
	void SimObserver::Update(const SimEvent& e)
	{
		// Only lock and proceed if we have to.
		if(e.GetIteration() % _frequency == 0 || e.ForceObserve())
		{
			_event = e;
			PreVisit();
			_event.GetObservable()->AcceptVisitor(*this);
			PostVisit();
		}
	};

	SimObserver* SimObserver::BuildObserver(const Json::Value &json,
									boost::mpi::communicator world,
									boost::mpi::communicator comm,
									int numwalkers,
									int wid)
	{
		return BuildObserver(json, world, comm, numwalkers, wid, "#/Observers");
	}

	SimObserver* SimObserver::BuildObserver(const Json::Value &json,
								boost::mpi::communicator world,
								boost::mpi::communicator comm,
								int numwalkers,
								int wid,
								const std::string& path)
	{
		ObjectRequirement validator; 
		Value schema;
		Reader reader;

		SimObserver* obs = nullptr;

		auto type = json.get("type", "none").asString();

		if(type == "JSON")
		{
			reader.parse(JsonSchema::JSONObserver, schema);
			validator.Parse(schema, path);

			// Validate inputs. 
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			// Initialize JSON specific options.
			auto prefix = json.get("file name", "unnamed").asString();
			auto frequency = json.get("frequency", 1).asInt();

			auto* dlm = new JSONObserver(prefix, world, comm, numwalkers, wid, frequency);
			obs = static_cast<SimObserver*>(dlm);
		}
		else
		{
			throw BuildException({path + ": Unknown observer type specified."});
		}

		return obs;
	}

	void SimObserver::BuildObservers(const Json::Value &json,
								boost::mpi::communicator world,
								boost::mpi::communicator comm,
								int numwalkers,
								int wid,
								ObserverList &ol)
	{
		ArrayRequirement validator;
		Value schema;
		Reader reader;

		reader.parse(JsonSchema::Observers, schema);
		validator.Parse(schema, "#/Observers");

		// Validate high level schema.
		validator.Validate(json, "#/Observers");
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		// Loop through observers.
		int i = 0;
		for(auto& obs : json)
		{
			ol.push_back(BuildObserver(obs, world, comm, numwalkers, wid, "#/Observers/" + std::to_string(i)));
			++i;
		}
	}
}
