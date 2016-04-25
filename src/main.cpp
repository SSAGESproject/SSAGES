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

	Simulation Sim(world);

	// Perform all the JSON reading and build the correct driver
	Sim.BuildSimulation(argv[1]);
	Sim.BuildDriver();

	// Run the MDEngine with Free energy calculations :)
	Sim.Run();
}
