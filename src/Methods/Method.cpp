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
#include "../Validator/ObjectRequirement.h"
#include "../Validator/ArrayRequirement.h"
#include "ABF.h"
#include "Umbrella.h"
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

		if(type == "Umbrella")
		{
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

			auto name = json.get("file name","none").asString();

			Umbrella* m = nullptr;
			if(timesteps == 0)
				m = new Umbrella(world, comm, ksprings, centers0, name, freq);
			else
				m = new Umbrella(world, comm, ksprings, centers0, centers1, timesteps, name, freq);

			if(json.isMember("iteration"))
				m->SetIteration(json.get("iteration",0).asInt());

			if(json.isMember("log every"))
				m->SetLogStep(json.get("log every",0).asInt());

			method = static_cast<Method*>(m);
		}
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
		else if(type == "ABF")
		{
			reader.parse(JsonSchema::ABFMethod, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<double> minsCV;
			for(auto& mins : json["CV_lower_bounds"])
				minsCV.push_back(mins.asDouble());
			
			std::vector<double> maxsCV;
			for(auto& maxs : json["CV_upper_bounds"])
				maxsCV.push_back(maxs.asDouble());

			std::vector<double> binsCV;
			for(auto& bins : json["CV_bins"])
				binsCV.push_back(bins.asDouble());

			std::vector<double> minsrestCV;
			for(auto& mins : json["CV_restraint_minimums"])
				minsrestCV.push_back(mins.asDouble());
			
			std::vector<double> maxsrestCV;
			for(auto& maxs : json["CV_restraint_maximums"])
				maxsrestCV.push_back(maxs.asDouble());

			std::vector<double> springkrestCV;
			for(auto& bins : json["CV_restraint_spring_constants"])
				springkrestCV.push_back(bins.asDouble());

			std::vector<bool> isperiodic;
			for(auto& isperCV : json["CV_isperiodic"])
				isperiodic.push_back(isperCV.asBool());

			std::vector<double> minperboundaryCV;
			for(auto& minsperCV : json["CV_periodic_boundary_lower_bounds"])
				minperboundaryCV.push_back(minsperCV.asDouble());

			std::vector<double> maxperboundaryCV;
			for(auto& maxsperCV : json["CV_periodic_boundary_upper_bounds"])
				maxperboundaryCV.push_back(maxsperCV.asDouble());

			
			if(!(minsCV.size() 	 	== maxsCV.size() && 
			     maxsCV.size() 	 	== binsCV.size() &&
			     binsCV.size()	 	== minsrestCV.size() &&
			     minsrestCV.size()	 	== maxsrestCV.size() &&
			     maxsrestCV.size()    	== springkrestCV.size() &&
			     springkrestCV.size()	== isperiodic.size()))			

			throw BuildException({"CV lower bounds, upper bounds, bins, restrain minimums, restrains maximums, spring constants and periodicity info must all have the size == number of CVs defined."});

			bool anyperiodic=false;
			for(size_t i = 0; i<isperiodic.size(); ++i)
				if(isperiodic[i])
					{
					anyperiodic = true;
					if(!(isperiodic.size() 	    	== minperboundaryCV.size() &&
			     		   minperboundaryCV.size()	== maxperboundaryCV.size()))
					throw BuildException({"If any CV is defined as periodic, please define the full upper and lower bound vectors. They should both have the same number of entries as CV lower bounds, upper bounds... Entries corresponding to non-periodic CVs will not be used."});
					}
			     
			int FBackupInterv = json.get("backup_frequency", 1000).asInt();

			double unitconv = json.get("unit_conversion", 0).asDouble();
		
			double timestep = json.get("timestep",2).asDouble();

			double min = json.get("minimum_count",100).asDouble();

			bool massweigh = json.get("mass_weighing",false).asBool();			

			std::vector<std::vector<double>> histdetails;
			std::vector<std::vector<double>> restraint;
			std::vector<std::vector<double>> periodicboundaries;
			std::vector<double> temp1(3);
			std::vector<double> temp2(3);
			std::vector<double> temp3(2);

			for(size_t i=0; i<minsCV.size(); ++i)
				{
				temp1 = {minsCV[i], maxsCV[i], binsCV[i]};
				temp2 = {minsrestCV[i], maxsrestCV[i], springkrestCV[i]};
				histdetails.push_back(temp1);
				restraint.push_back(temp2);
				if(anyperiodic)
					{
					temp3 = {minperboundaryCV[i], maxperboundaryCV[i]};
					periodicboundaries.push_back(temp3);
					}
				}

			auto freq = json.get("frequency", 1).asInt();

			std::string filename = json.get("filename", "F_out").asString();

			auto* m = new ABF(world, comm, restraint, isperiodic, periodicboundaries, min, massweigh, filename, histdetails, FBackupInterv, unitconv, timestep, freq);

			method = static_cast<Method*>(m);

			if(json.isMember("F") && json.isMember("N"))
				{
				Eigen::VectorXd F;
				std::vector<int> N;
				F.resize(json["F"].size());
				for(int i = 0; i < (int)json["F"].size(); ++i)
					F[i] = json["F"][i].asDouble();

				for(auto& n : json["N"])
					N.push_back(n.asInt());

				m->SetHistogram(F,N);
				}

			if(json.isMember("iteration"))
				m->SetIteration(json.get("iteration",0).asInt());

		}
		else if(type == "ForwardFlux")
		{
			reader.parse(JsonSchema::ForwardFluxMethod, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
            

            double ninterfaces = json.get("nInterfaces", 2).asDouble();
            
            std::vector<double> interfaces;
            for(auto& s : json["interfaces"])
            	interfaces.push_back(s.asDouble());

            std::vector<unsigned int> M;
            for(auto& s : json["trials"])
            	M.push_back(s.asInt());
           
           	if ((ninterfaces != interfaces.size()) || (ninterfaces != M.size())){
           		throw BuildException({"The size of \"interfaces\" and \"trials\" must be equal to \"nInterfaces\". See documentation for more information"});
           	}

            unsigned int N0Target = json.get("N0Target", 1).asInt();

            bool initialFluxFlag = json.get("computeInitialFlux", true).asBool();
           
            bool saveTrajectories = json.get("saveTrajectories", true).asBool();

            unsigned int currentInterface = json.get("currentInterface", 0).asInt();

            unsigned int freq = json.get("frequency", 1).asInt();

            std::string flavor = json.get("flavor", "none").asString();

            std::string output_directory = json.get("outputDirectoryName", "FFSoutput").asString();

            if(flavor == "DirectForwardFlux"){
            	auto *m = new DirectForwardFlux(world, comm, ninterfaces, interfaces, N0Target, M, initialFluxFlag, saveTrajectories, currentInterface, output_directory, freq);            	
            	method = static_cast<Method*>(m);
            } else {
            	throw BuildException({"Unknow flavor of forward flux. The options are \"DirectForwardFlux\""});
            }

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

