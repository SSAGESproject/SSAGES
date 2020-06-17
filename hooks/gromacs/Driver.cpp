/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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

#include "Driver.h"
#include "GromacsHook.h"
#include "ResourceHandler.h"
#include "json/json.h"
#include "gmxpre.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "programs/mdrun/mdrun_main.h"
#include <sstream>

// Starting in version 2019, gmx_mdrun is now in the gmx namespace.
// This allows us to call either &gmx_mdrun or &gmx::gmx_mdrun, depending on version.
using namespace gmx;

namespace SSAGES
{
	void Driver::Run()
	{
		int argc = args_.size() + 1 + ((rh_->GetNumWalkers() > 1) ? 2 : 0);
		char **argv = new char*[argc];
		for(int i = 0; i < argc; ++i)
			argv[i] = new char[128];
		
		sprintf(argv[0], "ssages");
		for(size_t i = 0; i < args_.size(); ++i)
			strcpy(argv[i+1], args_[i].c_str());
		
		if(rh_->GetNumWalkers() > 1)
		{
			sprintf(argv[args_.size()+1], "-multi");
			sprintf(argv[args_.size()+2], "%i", static_cast<int>(rh_->GetNumWalkers()));
		}

		std::cout << std::endl;
		CommandLineModuleManager::runAsMainCMain(argc, argv, &gmx_mdrun);
	}
	
	Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
	{
		auto* rh = ResourceHandler::Build(json, world);
		auto& hook = GromacsHook::Instance();
		rh->ConfigureHook(dynamic_cast<Hook*>(&hook));

		std::vector<std::string> args;
		if(json["args"].isArray())
		{
			for(auto& s : json["args"])
				args.push_back(s.asString());
		}
		else
		{
			std::istringstream iss(json["args"].asString());
			args = std::vector<std::string>(std::istream_iterator<std::string>{iss},
			                                std::istream_iterator<std::string>());
		}

		return new Driver(rh, args);
	}

	Driver::~Driver()
	{
		delete rh_;
	}
}
