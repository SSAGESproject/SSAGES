/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 Copyright 2017 Ben Sikora <bsikora906@gmail.com>
 *              Ashley Guo <azguo@uchicago.edu>
 *              Cody Bezik <bezik@uchicago.edu>
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

#include "ElasticBand.h"
#include "FiniteTempString.h"
#include "StringMethod.h"
#include "Swarm.h"
#include "CVs/CollectiveVariable.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "spline.h"
#include "schema.h"

using namespace Json;

namespace SSAGES
{
    void StringMethod::PrintString(const CVList& CV)
    {
        if(comm_.rank() == 0)
        {
            //Write node, iteration, centers of the string and current CV value to output file
            stringout_.precision(8);
            stringout_ << mpiid_ << " " << iteration_ << " ";

            for(size_t i = 0; i < centers_.size(); i++)
                stringout_ << worldstring_[mpiid_][i] << " " << CV[i]->GetValue() << " ";

            stringout_ << std::endl;
        }
    }

    void StringMethod::GatherNeighbors(std::vector<double> *lcv0, std::vector<double> *ucv0)
    {
        MPI_Status status;

        if(comm_.rank() == 0)
        {
            MPI_Sendrecv(&centers_[0], centers_.size(), MPI_DOUBLE, sendneigh_, 1234,
                &(*lcv0)[0], centers_.size(), MPI_DOUBLE, recneigh_, 1234, 
                world_, &status);

            MPI_Sendrecv(&centers_[0], centers_.size(), MPI_DOUBLE, recneigh_, 4321,
                &(*ucv0)[0], centers_.size(), MPI_DOUBLE, sendneigh_, 4321, 
                world_, &status);
        }

        MPI_Bcast(&(*lcv0)[0],centers_.size(),MPI_DOUBLE,0,comm_);
        MPI_Bcast(&(*ucv0)[0],centers_.size(),MPI_DOUBLE,0,comm_);
    }

    void StringMethod::StringReparam(double alpha_star)
    {
        std::vector<double> alpha_star_vector(numnodes_,0.0);

        //Reparameterization
        //Alpha star is the uneven mesh, approximated as linear distance between points
        if(comm_.rank()==0)
            alpha_star_vector[mpiid_] = mpiid_ == 0 ? 0 : alpha_star;

        //Gather each alpha_star into a vector 
        MPI_Allreduce(MPI_IN_PLACE, &alpha_star_vector[0], numnodes_, MPI_DOUBLE, MPI_SUM, world_);

        for(size_t i = 1; i < alpha_star_vector.size(); i++)
            alpha_star_vector[i] += alpha_star_vector[i-1];
        
        for(size_t i = 1; i < alpha_star_vector.size(); i++)
            alpha_star_vector[i] /= alpha_star_vector[numnodes_ - 1];

        tk::spline spl; //Cubic spline interpolation

        for(size_t i = 0; i < centers_.size(); i++)
        {
            std::vector<double> cvs_new(numnodes_, 0.0);

            if(comm_.rank() == 0)
                cvs_new[mpiid_] = centers_[i];

            MPI_Allreduce(MPI_IN_PLACE, &cvs_new[0], numnodes_, MPI_DOUBLE, MPI_SUM, world_);

            spl.set_points(alpha_star_vector, cvs_new);
            centers_[i] = spl(mpiid_/(numnodes_ - 1.0)); 
        }
    }

    void StringMethod::UpdateWorldString(const CVList& cvs)
    {
        for(size_t i = 0; i < centers_.size(); i++)
        {
            std::vector<double> cvs_new(numnodes_, 0.0);

            if(comm_.rank() == 0)
            {
                cvs_new[mpiid_] = centers_[i];
            }

            MPI_Allreduce(MPI_IN_PLACE, &cvs_new[0], numnodes_, MPI_DOUBLE, MPI_SUM, world_);

            for(int j = 0; j < numnodes_; j++)
            {
                worldstring_[j][i] = cvs_new[j];
                //Represent worldstring in periodic space
                worldstring_[j][i] = cvs[i]->GetPeriodicValue(worldstring_[j][i]); 
            }
        }
    }
    
    bool StringMethod::CheckEnd(const CVList& CV) 
    {
        if(maxiterator_ && iteration_ > maxiterator_)
        {
            std::cout << "System has reached max string method iterations (" << maxiterator_ << ") as specified in the input file(s)." << std::endl; 
            std::cout << "Exiting now" << std::endl; 
            PrintString(CV); //Ensure that the system prints out if it's about to exit
            MPI_Abort(world_, EXIT_FAILURE);
        }

        int local_tolvalue = TolCheck();

        MPI_Allreduce(MPI_IN_PLACE, &local_tolvalue, 1, MPI_INT, MPI_LAND, world_);

        if(local_tolvalue)
        {
            std::cout << "System has converged within tolerance criteria. Exiting now" << std::endl;
            PrintString(CV); //Ensure that the system prints out if it's about to exit
            MPI_Abort(world_, EXIT_FAILURE);
        }

        return true;
    }

    void StringMethod::PreSimulation(Snapshot* snapshot, const CVList& cvs)
    {
        mpiid_ = snapshot->GetWalkerID();
        char file[1024];
        sprintf(file, "node-%04d.log",mpiid_);
        stringout_.open(file);

        SetSendRecvNeighbors();

        worldstring_.resize(numnodes_);
        for(auto& w : worldstring_)
            w.resize(centers_.size());
        UpdateWorldString(cvs);
        PrintString(cvs);
    }

    void StringMethod::SetSendRecvNeighbors()
    {
        std::vector<int> wiids(world_.size(), 0);

        //Set the neighbors
        recneigh_ = -1;
        sendneigh_ = -1; 

        MPI_Allgather(&mpiid_, 1, MPI_INT, &wiids[0], 1, MPI_INT, world_);
        numnodes_ = int(*std::max_element(wiids.begin(), wiids.end())) + 1;

        // Ugly for now...
        for(size_t i = 0; i < wiids.size(); i++)
        {
            if(mpiid_ == 0)
            {
                sendneigh_ = comm_.size();
                if(wiids[i] == numnodes_ - 1)
                {
                    recneigh_ = i;
                    break;
                }
            }
            else if (mpiid_ == numnodes_ - 1)
            {
                sendneigh_ = 0;
                if(wiids[i] == mpiid_ - 1)
                {
                    recneigh_ = i;
                    break;
                }
            } 
            else
            {
                if(wiids[i] == mpiid_ + 1)
                {
                    sendneigh_ = i;
                    break;
                }
                if(wiids[i] == mpiid_ - 1 && recneigh_ == -1)
                    recneigh_ = i;
            }
        }
    }

    void StringMethod::Serialize(Json::Value& json) const
    {
        json["type"] = "String";

        for(size_t i = 0; i < centers_.size(); i++)
            json["centers"].append(centers_[i]);

        for(auto& t : tol_)
            json["tolerance"].append(t);

        json["max_iterations"] = maxiterator_;

        for(auto& s : cvspring_)
            json["ksprings"].append(s);

        json["iteration"] = iteration_; 
    }

	//! \copydoc Buildable::Build()
	StringMethod* StringMethod::Construct(const Value& json, 
		                                  const MPI_Comm& world,
		                                  const MPI_Comm& comm,
					                      const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;

		StringMethod* m = nullptr;

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
		if(flavor == "ElasticBand")
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

			m = new ElasticBand(world, comm, centers, 
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
		}
		else if(flavor == "FTS")
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
			m = new FiniteTempString(world, comm, centers, 
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
			
			if(json.isMember("tolerance"))
			{
				std::vector<double> tol;
				for(auto& s : json["tolerance"])
					tol.push_back(s.asDouble());

				m->SetTolerance(tol);
			}
		}
		
		return m;
	}

}