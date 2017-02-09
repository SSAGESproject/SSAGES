/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Hythem Sidky <hsidky@nd.edu>
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

#include "lammps.h"
#include "../Observers/JSONObserver.h"
#include "Drivers/Driver.h"
#include "Drivers/DriverException.h"
#include "../Validator/ObjectRequirement.h"
#include "schema.h"
#include "input.h"
#include "modify.h"
#include "output.h"
#include "update.h"
#include "fix.h"

namespace SSAGES
{
	//! Driver for LAMMPS simulations
	class LammpsDriver : public Driver 
	{
	private:

		//! Pointer to the local instance of lammps
		std::shared_ptr<LAMMPS_NS::LAMMPS> lammps_;

		//! The lammps logfile
		std::string logfile_;

	public:

		//! Constructor
		/*!
		 * \param world_comm MPI global communicator.
		 * \param local_comm MPI local communicator.
		 * \param walkerID ID of the walker assigned to this driver.
		 */
		LammpsDriver(boost::mpi::communicator& world_comm,
					 boost::mpi::communicator& local_comm,
					 int walkerID) : 
		Driver(world_comm, local_comm, walkerID), lammps_(), logfile_() 
		{
		};

		//! Run simulation
		virtual void Run() override
		{
			std::string rline = "run " + std::to_string(iterations_);
			lammps_->input->one(rline.c_str());
		}

		//! Run LAMMPS input file
		/*!
		 * \param contents Content of the LAMMPS input file.
		 *
		 * This function exectures the contents of the given LAMMPS input file
		 * line by line and gathers the fix/hook.
		 */
		virtual void ExecuteInputFile(std::string contents) override
		{
			// Go through lammps.
			std::uniform_int_distribution<> dis(1, 999999);

			std::string token;
			std::istringstream ss(contents);
			bool reading_restart = false;

			if(restartname_ != "none" && readrestart_)
			{
				lammps_->input->one(("read_restart " + restartname_ + " remap").c_str());
				while(std::getline(ss, token, '\n'))
					if(token.find("#RESTART") != std::string::npos)
						reading_restart = true;

			}

			// need these 2 lines
			ss.clear(); // clear the `failbit` and `eofbit`
			ss.seekg(0); // rewind

			while(std::getline(ss, token, '\n'))
			{
				if(token.find("#RESTART") != std::string::npos)
					reading_restart = false;

				if(readrestart_ && reading_restart)
					continue;

				lammps_->input->one(token.c_str());
			}

			// Initialize and create the restart parameters
			for(auto* o : observers_)
			{
				if(o->GetName() == "JSON")
				{
					readrestart_ = true;
					JSONObserver* obs = static_cast<JSONObserver*>(o);
					std::string filename1 = obs->GetPrefix() + "_" + std::to_string(wid_) + ".restart";
					std::string filename2 = obs->GetPrefix() + "_" + std::to_string(wid_) + "b.restart";
					lammps_->input->one(("restart  " + std::to_string(obs->GetFrequency()) + " " + filename1 + " " + filename2).c_str());
				}
			}

			auto fid = lammps_->modify->find_fix("ssages");
			if(fid < 0)
				throw BuildException({"Could not find ssages fix in given input file!"});

			if(!(hook_ = dynamic_cast<Hook*>(lammps_->modify->fix[fid])))
			{
				throw BuildException({"Unable to dynamic cast hook on node " + std::to_string(world_.rank())});			
			}
		}
		
		//! Set up the driver
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		virtual void BuildDriver(const Json::Value& json, const std::string& path) override
		{

			Json::Value schema;
			Json::ObjectRequirement validator;
			Json::Reader reader;

			reader.parse(JsonSchema::LAMMPSDriver, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			iterations_ = json.get("MDSteps",1).asInt();
			inputfile_ = json.get("inputfile","none").asString();
			restartname_ = json.get("restart file", "none").asString();
			readrestart_ = json.get("read restart", false).asBool();

			if(readrestart_ && restartname_ == "none")
				throw BuildException({"You want to run from a restart but no file name provided (see 'restart file' in LAMMPS's schema for more informationz)"});
			
			// Silence of the lammps.
			char **largs = (char**) malloc(sizeof(char*) * 1);
			largs[0] = (char*) malloc(sizeof(char) * 1024);
			sprintf(largs[0], " ");
			lammps_ = std::make_shared<LAMMPS_NS::LAMMPS>(1, largs, MPI_Comm(comm_));

			// Free.
			for(int i = 0; i < 1; ++i)
				free(largs[i]);
			free(largs);
		}

		virtual void Serialize(Json::Value& json) const override
		{
			// Call parent first.
			Driver::Serialize(json);

			json["MDSteps"] = iterations_;
			json["logfile"] = logfile_;
			json["type"] = "LAMMPS";
			//if true on first file
			if(lammps_->output->restart_toggle)
			{
				json["restart file"] = std::string(lammps_->output->restart2a);
			}
			else
			{
				json["restart file"] = std::string(lammps_->output->restart2b);
			}
		}
	};
}
