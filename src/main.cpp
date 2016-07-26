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

	if(world.rank()==0)
		std::cout << std::setw(47 + 8) << std::left << "\033[1m > JSON validation and simulation building\n \033[0m" << std::flush;

	bool SimCheck = Sim.BuildSimulation(root, "#/Simulations");
	bool DriverCheck = Sim.BuildDriver(root, "#/Drivers");

	if(DriverCheck == false || SimCheck == false)
	{
		std::cout<<"Simulation/Grid/Method/CV/Driver build fail"<<std::endl;
		world.abort(-1);
	}
	world.barrier();
	
	try{
		if(world.rank()==0)
			std::cout << std::setw(47 + 8) << std::left << "\033[1m > Reading MD Engine Input file...\033[0m" << std::endl;	
		Sim.ReadInputFile();
		if(world.rank()==0)
			std::cout << std::setw(34) << std::right << "\033[32mFinished!\033[0m\n";
		world.barrier();
		
		if(world.rank()==0)
			std::cout << std::setw(47 + 8) << std::left << "\033[1m > Finalizing Objects\033[0m" << std::flush;
		Sim.Finalize();
		if(world.rank()==0)
			std::cout << std::setw(34) << std::right << "\033[32mFinished!\033[0m\n";
		world.barrier();
		
		// Run the MDEngine with Free energy calculations :)
		if(world.rank()==0)
			std::cout << std::setw(47 + 8) << std::left << "\033[1m > Running simulation... \033[0m" << std::flush;
		Sim.Run();
		if(world.rank()==0)
			std::cout << std::setw(34) << std::right << "\033[32mFinished!\033[0m\n";
	} catch(BuildException& e) {
		DumpErrorsToConsole(e.GetErrors(), 30);
		world.abort(-1);
	} catch(std::exception& e) {
		DumpErrorsToConsole({e.what()}, 30);
		world.abort(-1);
	} catch(int& k) { 
		std::string err = strerror(k);
		DumpErrorsToConsole({"File IO error: " + err}, 30);
		world.abort(-1);
	}
}
