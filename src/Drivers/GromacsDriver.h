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
	public:
		GromacsDriver(mpi::communicator& world_comm, 
					  mpi::communicator& local_comm,
					  int walkerID) : 
		Driver(world_comm, local_comm, walkerID),
		nwalkers_(world_comm.size()/local_comm.size())
		{}

		virtual void Run() override
		{
			int argc = 5; 
			char **largs = new char*[argc];
			largs[0] = new char[128];
			largs[1] = new char[128];
			largs[2] = new char[128];
			largs[3] = new char[128];
			largs[4] = new char[128];

			// Trim input file extension.
			auto s = _inputfile.substr(0, _inputfile.find_last_of("."));

			sprintf(largs[0], "ssages");
			sprintf(largs[1], "-multi");
			sprintf(largs[2], "%i", nwalkers_);
			sprintf(largs[3], "-deffnm");
			sprintf(largs[4], "%s", s.c_str());

			// For prettyness.
			std::cout << std::endl;
			gmx::CommandLineModuleManager::runAsMainCMain(argc, largs, &gmx_mdrun);
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
			auto iterations = json.get("MDSteps", 1).asInt();

			// Set hook. 
			auto& hook = GromacsHook::Instance();
			hook.SetIterationTarget(iterations);
			_hook = dynamic_cast<Hook*>(&hook);
		}

		virtual void Serialize(Json::Value& json) const override
		{

		}
	};
}