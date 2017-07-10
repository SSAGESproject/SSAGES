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

#include "Driver.h"
#include "GromacsHook.h"
#include "ResourceHandler.h"
#include "gmxpre.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "programs/mdrun/mdrun_main.h"

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
			sprintf(argv[args_.size()+2], "%i", (int)rh_->GetNumWalkers());
		}

		std::cout << std::endl;
		gmx::CommandLineModuleManager::runAsMainCMain(argc, argv, &gmx_mdrun);
	}
	
	Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
	{
		auto* rh = ResourceHandler::Build(json, world);
		auto& hook = GromacsHook::Instance();
		rh->ConfigureHook(dynamic_cast<Hook*>(&hook));

		std::vector<std::string> args; 
		for(auto& s : json["args"])
			args.push_back(s.asString());

		return new Driver(rh, args);
	}   

	Driver::~Driver()
	{
		delete rh_;
	} 
}
