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

	Simulation Sim(world);

	// Perform all the JSON reading and build the correct driver
	root = Sim.ReadJSON(argv[1]);
	Sim.BuildSimulation(root, "#/Simulations");
	Sim.BuildDriver(root.get("driver", Json::arrayValue), "#/Drivers");
	Sim.BuildCVs();
	Sim.BuildMethod();

	// Run the MDEngine with Free energy calculations :)
	Sim.Run();
}
