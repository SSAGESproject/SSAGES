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

		else if(type == "Metadynamics")
		{
			reader.parse(JsonSchema::MetadynamicsMethod, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<double> widths;
			for(auto& s : json["widths"])
				widths.push_back(s.asDouble());

			std::vector<double> lowerb, upperb, lowerk, upperk;
			Grid<Vector>* grid = nullptr;
			if(json.isMember("grid"))
				grid = Grid<Vector>::BuildGrid(json.get("grid", Json::Value()));
			else if(!json.isMember("lower_bounds") || !json.isMember("upper_bounds"))
				throw BuildException({
					"#/Method/Metadynamics: Both upper_bounds and lower_bounds "
					"must be defined if grid is not being used."});

			// Assume all vectors are the same size. 
			for(int i = 0; i < json["lower_bound_restraints"].size(); ++i)
			{
				lowerk.push_back(json["lower_bound_restraints"][i].asDouble());
				upperk.push_back(json["upper_bound_restraints"][i].asDouble());
				lowerb.push_back(json["lower_bounds"][i].asDouble());
				upperb.push_back(json["upper_bounds"][i].asDouble());
			}
		
			auto height = json.get("height", 1.0).asDouble();
			auto hillfreq = json.get("hill_frequency", 1).asInt();
			auto freq = json.get("frequency", 1).asInt();

			auto* m = new Meta(
			    world, comm, height, widths, 
				lowerb, upperb, lowerk,	upperk,
				grid, hillfreq, freq
			);

			if(json.isMember("load_hills"))
				m->LoadHills(json["load_hills"].asString());

			method = static_cast<Method*>(m);
		}
		
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
        else if(type == "Basis")
        {
            reader.parse(JsonSchema::BFSMethod, schema);
            validator.Parse(schema, path);

            //Validate Inputs
            validator.Validate(json, path);
            if(validator.HasErrors())
                throw BuildException(validator.GetErrors());

			std::vector<unsigned int> coefsCV(0);
			for(auto& coefs : json["CV_coefficients"])
				coefsCV.push_back(coefs.asInt());

            std::vector<double> restrCV(0);
			for(auto& restr : json["CV_restraint_spring_constants"])
                restrCV.push_back(restr.asDouble());

            std::vector<double> boundLow(0);
            for(auto& bndl : json["CV_restraint_minimums"])
                boundLow.push_back(bndl.asDouble());
        
            std::vector<double> boundUp(0);
            for(auto& bndu : json["CV_restraint_maximums"])
                boundUp.push_back(bndu.asDouble());

            auto cyclefreq = json.get("cycle_frequency", 100000).asInt();
            auto freq = json.get("frequency", 1).asInt();
            auto wght = json.get("weight", 1.0).asDouble();
            auto bnme = json.get("basis_filename", "").asString();
            auto cnme = json.get("coeff_filename", "").asString();
            auto temp = json.get("temperature", 0.0).asDouble();
            auto tol  = json.get("tolerance", 1e-6).asDouble();
            auto conv = json.get("convergence_exit", false).asBool();
 
            Histogram<int> *hist = Histogram<int>::BuildHistogram(
                                            json.get("grid", Json::Value()) );

            auto* m = new Basis(world, comm, hist, coefsCV, restrCV, boundUp, boundLow,
                                cyclefreq, freq, bnme, cnme, temp, tol, wght,
                                conv);

            method = static_cast<Method*>(m);
			
            if(json.isMember("iteration"))
				m->SetIteration(json.get("iteration",0).asInt());

            if(json.isMember("coefficients") && json.isMember("bias hist"))
            {
                std::vector<double> coeff;
                std::vector<double> unbias;

                for(auto& c : json["coefficients"])
                    coeff.push_back(c.asDouble());

                for(auto& u : json["bias hist"])
                    unbias.push_back(u.asDouble());

                m->SetBasis(coeff, unbias);
            }

        }
		else
		{
			throw BuildException({path + ": Unknown method type specified."});
		}
		
		return method;
	}
	*/
}

