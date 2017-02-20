/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
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

#include "gmxpre.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "programs/mdrun/mdrun_main.h"
#include "Driver.h"
#include "../../hooks/gromacs/GromacsHook.h"

namespace SSAGES
{
	//! Driver for Gromacs simulations.
	class GromacsDriver : public Driver
	{
	private:

		//! Number of walkers.
		int nwalkers_;

		//! Is the simulation a restart?
		bool restart_;

	public:

		//! Constructor.
		/*!
		 * \param world_comm MPI global communicator.
		 * \param local_comm MPI local communicator.
		 * \param walkerID ID of the walker assigned to this driver.
		 */
		GromacsDriver(mpi::communicator& world_comm, 
					  mpi::communicator& local_comm,
					  int walkerID) : 
		Driver(world_comm, local_comm, walkerID),
		nwalkers_(world_comm.size()/local_comm.size())
		{}

		//! Run simulation.
		virtual void Run() override
		{
			int argc = 3; 
			char **largs = new char*[64];
			largs[0] = new char[128];
			largs[1] = new char[128];
			largs[2] = new char[128];

			// Trim input file extension.
			auto s = inputfile_.substr(0, inputfile_.find_last_of("."));

			sprintf(largs[0], "ssages");
			sprintf(largs[1], "-deffnm");
			sprintf(largs[2], "%s", s.c_str());

			if(nwalkers_ > 1)
			{
				largs[argc] = new char[128];
				sprintf(largs[argc], "-multi");
				++argc;
				largs[argc] = new char[128];
				sprintf(largs[argc], "%i", nwalkers_);
				++argc;

			}

			if(restart_)
			{
				largs[argc] = new char[128];
				sprintf(largs[argc], "-cpi");
				++argc;
				largs[argc] = new char[128];
				sprintf(largs[argc], "-noappend");
				++argc;
			}

			// For prettyness.
			std::cout << std::endl;
			gmx::CommandLineModuleManager::runAsMainCMain(argc, largs, &gmx_mdrun);
		}

		//! Run Input file.
		/*!
		 * In the case of Gromacs, the input file does not need to be parsed for
		 * special information.
		 */
		virtual void ExecuteInputFile(std::string) override
		{
		}

		//! Set up the driver.
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		virtual void BuildDriver(const Json::Value& json, const std::string& path) override
		{
			inputfile_ = json.get("inputfile","none").asString();
			iterations_ = json.get("MDSteps", 1).asInt();

			// Set hook. 
			auto& hook = GromacsHook::Instance();
			hook.SetIterationTarget(iterations_);
			hook_ = dynamic_cast<Hook*>(&hook);

			// Restart?
			restart_ = json.get("restart", false).asBool();
		}

		//! \copydoc Serializable::Serialize()
		virtual void Serialize(Json::Value& json) const override
		{
			// Call parent first.
			Driver::Serialize(json);

			json["MDSteps"] = iterations_;
			json["type"] = "Gromacs";
		}
	};
}
