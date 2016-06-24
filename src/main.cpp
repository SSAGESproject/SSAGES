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
	bool SimCheck = Sim.BuildSimulation(root, "#/Simulations");
	bool CVCheck = Sim.BuildCVs(root, "#/CVs");
	bool MethodCheck = Sim.BuildMethod(root, "#/Methods");
	bool DriverCheck = Sim.BuildDriver(root.get("driver", Json::arrayValue), "#/Drivers");
	bool GridCheck = Sim.BuildGrid(root, "#/Grids");

	if(CVCheck == false || MethodCheck == false || DriverCheck == false || SimCheck == false || GridCheck == false)
	{
		std::cout<<"Simulation/Grid/Method/CV/Driver build fail"<<std::endl;
		world.abort(-1);
	}

	std::cout << std::setw(47 + 8) << std::left << "\033[1m > JSON validation and simulation building \033[0m" << std::flush;
	std::cout << std::setw(34) << std::right << "\033[32mFinished!\033[0m\n";
	std::cout << std::setw(47 + 8) << std::left << "\033[1m > Reading MD Engine Input file...\033[0m" << std::endl;;	
	Sim.ReadInputFile();
	std::cout << std::setw(34) << std::right << "\033[32mFinished!\033[0m\n";

	std::cout << std::setw(47 + 8) << std::left << "\033[1m > Finalizing Objects\033[0m" << std::flush;
	Sim.Finalize();
	std::cout << std::setw(34) << std::right << "\033[32mFinished!\033[0m\n";

	// Run the MDEngine with Free energy calculations :)
	std::cout << std::setw(47 + 8) << std::left << "\033[1m > Running simulation... \033[0m" << std::flush;
	Sim.Run();
	std::cout << std::setw(34) << std::right << "\033[32mFinished!\033[0m\n";
}
