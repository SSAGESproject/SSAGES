#include <boost/mpi.hpp>
#include <string>
#include <sstream>
#include "json/json.h"
#include "JSON/JSONLoader.h"
#include "input.h"
#include "modify.h"
#include "Snapshot.h"
#include "Hook.h"
#include "fix.h"
#include "Validator/ObjectRequirement.h"
#include "Validator/ArrayRequirement.h"
#include "schema.h"
#include "Drivers/Driver.h"
#include "Utility.h"

namespace mpi = boost::mpi;

using namespace SSAGES;
using namespace Json;

int main(int argc, char* argv[])
{
	mpi::environment env(argc, argv);
	mpi::communicator world;

	Value root;
	ObjectRequirement validator;
	Reader reader;
	Driver* MDdriver = nullptr;

	// Read in JSON using head node and broadcast to other nodes
	// JSON file will include the Engine input file name(s)
	if(world.rank() == 0)
		root = ReadJSONFile(argv[1]);
	mpi::broadcast(world, root, 0);

	// Get requested number of walkers and make 
	// sure the number of processors is evenly divisible.
	int nwalkers = root.get("number walkers", 1).asInt();
	if(world.size() % nwalkers != 0)
	{
		if(world.rank() == 0)
			std::cerr << "The number of processors must be evenly "
			<< "divisible by the number of walkers." << std::endl;
		world.abort(-1);
	}

	std::string path = "#/Drivers";

	// Determine walker ID using integer division and split 
	// communicator accordingly. 
	int wid = (int)world.rank() / nwalkers;
	auto walker = world.split(wid);

	reader.parse(JsonSchema::Driver, root);
	validator.Parse(root, path);

	// Validate inputs.
	validator.Validate(root, path);
	if(validator.HasErrors())
		throw BuildException(validator.GetErrors());

	// Get move type. 
	std::string MDEngine = root.get("MDEngine", "none").asString();
	std::string inputfile = root.get("inputfile", "none").asString();

	// Use input from JSON to determine MDEngine of choice as well as other parameters
	if(MDEngine == "LAMMPS")
	{
		auto* en = new lammpsMD(world, walker, wid, root, inputfile);
		MDDriver = static_cast<Driver*>(en);
	}
	else
	{
		throw BuildException({"Unknown MD Engine specified."});
	}

	if(MDDriver == nullptr)
	{
		std::cerr << "Could not create MDDriver object!" << std::endl;
		world.abort(-1);
	}

	MDDriver->BuildDriver(root, "#/Drivers");

	// Read inputfile contents on rank 0 proc, broadcast to all.
	std::string contents;
	// All nodes get the same input file
	if( inputfile != "none")
	{
		if(world.rank() == 0)
			contents = GetFileContents(inputfile.c_str());
		mpi::broadcast(world, contents, 0);
	}
	// Each node reads in specified input file
	else(inputfile == "none")
	{
		if(walker.rank() == 0)
		{
			cout<<"No/overloaded global input file, node " <<wid<<" using: "<<MDEngine->GetInputFile()<<std::endl;
			contents = GetFileContents(MDEngine->GetInputFile())
		}
		mpi::broadcast(walker, contents, 0);
	}


	MDDriver->BuildCVs();
	MDDriver->BuildMethod();

	// Execute Global/local input file and gather the needed hook
	MDDriver->ExecuteInputFile(contents);

	// Setting up listeners
	MDDriver->Finalize();

	// Run the MDEngine with Free energy calculations :)
	MDDriver->Run();

	delete MDDriver;

}
