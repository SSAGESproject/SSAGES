#include <boost/mpi.hpp>
#include <string>
#include <sstream>
#include "Simulations/Simulation.h"


namespace mpi = boost::mpi;

using namespace SSAGES;
using namespace Json;

int main(int argc, char* argv[])
{
	mpi::environment env(argc, argv);
	mpi::communicator world;

	Json::Value root;
	Json::Value JsonDriver;

	Simulation Sim(world);

	// Perform all the JSON reading and build the correct driver
	root = Sim.ReadJSON(argv[1]);
	Sim.BuildSimulation(root, "#/Simulations");
	JsonDriver = Sim.BuildDriver(root.get("driver", Json::arrayValue), "#/Drivers");
	bool CVCheck = Sim.BuildCVs(JsonDriver, "#/CVs");
	bool MethodCheck = Sim.BuildMethod(JsonDriver, "#/Methods");

	if(CVCheck == 0 || MethodCheck == 0)
	{
		if(_comm.rank() == 0)
			std::cout<<"Method and/or CV fail on node "<<_world.rank()<<std::endl;
		_world.abort(-1);
	}

	Sim.ReadInputFile();
	Sim.Finalize();

	// Run the MDEngine with Free energy calculations :)
	Sim.Run();
}
