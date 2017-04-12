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
#include "Method.h"
#include "json/json.h"
#include "schema.h"
#include "Validator/ObjectRequirement.h"
#include "Validator/ArrayRequirement.h"
#include "ABF.h"
#include "Umbrella.h"
#include "BasisFunc.h"
#include "ForwardFlux.h"
#include "Meta.h"
#include "StringMethod.h"
#include <stdexcept>

using namespace Json;

namespace SSAGES
{
	Method* Method::BuildMethod(const Json::Value& json, 
		                        const MPI_Comm& world, 
							    const MPI_Comm& comm, 
							    const std::string& path)
	{
		if(json["type"] == "ABF")
			return ABF::Build(json, world, comm, path);
		else if(json["type"] == "Basis")
			return Basis::Build(json, world, comm, path);
		else if(json["type"] == "ForwardFlux")
			return ForwardFlux::Build(json, world, comm, path);
		else if(json["type"] == "Metadynamics")
			return Meta::Build(json, world, comm, path);
		else if(json["type"] == "Umbrella")
			return Umbrella::Build(json, world, comm, path);
		else if(json["type"] == "String")
			return StringMethod::Build(json, world, comm, path);
		else
			throw std::invalid_argument(path + ": Unknown method type specified.");
	}
	/*
	Method* Method::BuildMethod(const Value &json, 
						boost::mpi::communicator& world, 
						boost::mpi::communicator& comm,
						const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;

		Method* method = nullptr;

		// Get method type. 
		std::string type = json.get("type", "none").asString();

		else if(type == "String")
		{
			reader.parse(JsonSchema::StringMethod, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<double> centers;
			for(auto& s : json["centers"])
				centers.push_back(s.asDouble());

			auto maxiterator = json.get("max_iterations", 0).asInt();

			std::vector<double> ksprings;
			for(auto& s : json["ksprings"])
				ksprings.push_back(s.asDouble());

			auto freq = json.get("frequency", 1).asInt();

			// Get stringmethod flavor. 
			std::string flavor = json.get("flavor", "none").asString();

			if(flavor == "FTS")
			{
				reader.parse(JsonSchema::FTSMethod, schema);
				validator.Parse(schema, path);

				// Validate inputs.
				validator.Validate(json, path);
				if(validator.HasErrors())
					throw BuildException(validator.GetErrors());
    			auto isteps = json.get("block_iterations", 2000).asInt();
		    	auto tau = json.get("time_step", 0.1).asDouble();
				auto kappa = json.get("kappa", 0.1).asDouble();
				auto springiter = json.get("umbrella_iterations",2000).asDouble();
				auto* m = new FiniteTempString(world, comm, centers, 
									maxiterator, isteps,
									tau, ksprings, kappa,
									springiter, freq);

				if(json.isMember("tolerance"))
				{
					std::vector<double> tol;
					for(auto& s : json["tolerance"])
						tol.push_back(s.asDouble());

					m->SetTolerance(tol);
				}

				method = static_cast<Method*>(m);
			}
			else if(flavor == "ElasticBand")
			{
				reader.parse(JsonSchema::ElasticBandMethod, schema);
				validator.Parse(schema, path);

				// Validate inputs.
				validator.Validate(json, path);
				if(validator.HasErrors())
					throw BuildException(validator.GetErrors());

				auto eqsteps = json.get("equilibration_steps", 20).asInt();
				auto evsteps = json.get("evolution_steps", 5).asInt();
				auto stringspring = json.get("kstring", 10.0).asDouble();
    			auto isteps = json.get("block_iterations", 100).asInt();
		    	auto tau = json.get("time_step", 0.1).asDouble();

				auto* m = new ElasticBand(world, comm, centers, 
									maxiterator, isteps,
									tau, ksprings, eqsteps,
									evsteps, stringspring, freq);

                if(json.isMember("tolerance"))
				{
					std::vector<double> tol;
					for(auto& s : json["tolerance"])
						tol.push_back(s.asDouble());

					m->SetTolerance(tol);
				}

				method = static_cast<Method*>(m);
			}
            else if(flavor == "SWARM")
            {
                reader.parse(JsonSchema::SwarmMethod, schema);
                validator.Parse(schema, path);
                
                //Validate input
                validator.Validate(json, path);
                if(validator.HasErrors())
                    throw BuildException(validator.GetErrors());
 
                auto InitialSteps = json.get("initial_steps", 2500).asInt();
                auto HarvestLength = json.get("harvest_length", 10).asInt();
                auto NumberTrajectories = json.get("number_of_trajectories", 250).asInt();
                auto SwarmLength = json.get("swarm_length", 20).asInt();
                
                auto* m = new Swarm(world, comm, centers, maxiterator, ksprings, freq, InitialSteps, HarvestLength, NumberTrajectories, SwarmLength); 
                method = static_cast<Method*>(m);
            
                if(json.isMember("tolerance"))
				{
					std::vector<double> tol;
					for(auto& s : json["tolerance"])
						tol.push_back(s.asDouble());

					m->SetTolerance(tol);
				}

				method = static_cast<Method*>(m);
			}
			else
			{
				throw BuildException({flavor + " is unknown string method type. Please specify correct flavor"});
			}
            
            if(json.isMember("iteration"))
                method->SetIteration(json.get("iteration",0).asInt());
            
		}
		else
		{
			throw BuildException({path + ": Unknown method type specified."});
		}
		
		return method;
	}
	*/
}

