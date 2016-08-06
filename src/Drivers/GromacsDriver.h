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
	public:
		GromacsDriver(mpi::communicator& world_comm, 
					  mpi::communicator& local_comm,
					  int walkerID) : 
		Driver(world_comm, local_comm, walkerID)
		{}

		virtual void Run() override
		{
			int argc = 3; 
			char **largs = new char*[argc];
			largs[0] = new char[100];
			largs[1] = new char[100];
			largs[2] = new char[100];

			// Trim input file extension.
			auto s = _inputfile.substr(0, _inputfile.find_last_of("."));

			sprintf(largs[0], "ssages");
			sprintf(largs[1], "-deffnm");
			sprintf(largs[2], "%s", s.c_str());

			// For prettyness.
			std::cout << std::endl;
			gmx::CommandLineModuleManager::runAsMainCMain(argc, largs, &gmx_mdrun);

			for(int i = 0; i < argc; ++i)
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

			// Set hook. 
			_hook = dynamic_cast<Hook*>(&GromacsHook::Instance());
		}

		virtual void Serialize(Json::Value& json) const override
		{

		}
	};
}