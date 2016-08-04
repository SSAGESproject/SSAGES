#pragma once 

#include "gmxpre.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "programs/mdrun/mdrun_main.h"
#include "Driver.h"

namespace SSAGES
{
	class GromacsDriver : public Driver
	{
	public:
		GromacsDriver(mpi::communicator& world_comm, 
					  mpi::communicator& local_comm,
					  int walkerID) : 
		Driver(world_comm, local_comm, walkerID)
		{}

		virtual void Run() override
		{
			int argc = 2; 
			char **largs = new char*[2];
			largs[0] = new char[100];
			largs[1] = new char[100];
			sprintf(largs[0], "ssages");
			sprintf(largs[1], "help");
			gmx::CommandLineModuleManager::runAsMainCMain(argc, largs, &gmx_mdrun);

			for(int i = 0; i < 2; ++i)
				free(largs[i]);
			free(largs);
		}

		virtual void ExecuteInputFile(std::string contents) override
		{

		}

		virtual void WriteRestartFile() const override
		{

		}

		virtual void BuildDriver(const Json::Value& json, const std::string& path) override
		{
			_inputfile = json.get("inputfile","none").asString();
		}

		virtual void Serialize(Json::Value& json) const override
		{

		}
	};
}