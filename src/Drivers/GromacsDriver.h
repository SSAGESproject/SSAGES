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
	class GromacsDriver : public Driver
	{
	private:
		int nwalkers_;
		int MDSteps_;

	public:
		GromacsDriver(mpi::communicator& world_comm, 
					  mpi::communicator& local_comm,
					  int walkerID) : 
		Driver(world_comm, local_comm, walkerID),
		nwalkers_(world_comm.size()/local_comm.size())
		{}

		virtual void Run() override
		{
			int argc = (nwalkers_ > 1) ? 5 : 3; 
			char **largs = new char*[argc];
			largs[0] = new char[128];
			largs[1] = new char[128];
			largs[2] = new char[128];
			if(nwalkers_ > 1)
			{
				largs[3] = new char[128];
				largs[4] = new char[128];
			}

			// Trim input file extension.
			auto s = _inputfile.substr(0, _inputfile.find_last_of("."));

			sprintf(largs[0], "ssages");
			sprintf(largs[1], "-deffnm");
			sprintf(largs[2], "%s", s.c_str());
			if(nwalkers_ > 1)
			{
				sprintf(largs[3], "-multi");
				sprintf(largs[4], "%i", nwalkers_);
			}
				
			// For prettyness.
			std::cout << std::endl;
			gmx::CommandLineModuleManager::runAsMainCMain(argc, largs, &gmx_mdrun);
		}

		virtual void ExecuteInputFile(std::string) override
		{

		}

		virtual void BuildDriver(const Json::Value& json, const std::string&) override
		{
			_inputfile = json.get("inputfile","none").asString();
			MDSteps_ = json.get("MDSteps", 1).asInt();

			// Set hook. 
			auto& hook = GromacsHook::Instance();
			hook.SetIterationTarget(MDSteps_);
			_hook = dynamic_cast<Hook*>(&hook);
		}

		virtual void Serialize(Json::Value& json) const override
		{
			// Call parent first.
			Driver::Serialize(json);

			json["MDSteps"] = MDSteps_;
			json["type"] = "Gromacs";
		}
	};
}