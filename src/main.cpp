#include <boost/mpi.hpp>
#include <string>
#include <sstream>
#include "lammps.h"
#include "input.h"
#include "modify.h"
#include "Snapshot.h"
#include "Hook.h"
#include "fix.h"
#include "Methods/MockMethod.h"

namespace mpi = boost::mpi;
using namespace LAMMPS_NS;
using namespace SSAGES;

// Retrieves the contents of a file and returns them
// in a string. Throws exception on failure.
std::string GetFileContents(const char *filename)
{
	std::FILE *fp = std::fopen(filename, "rb");
	if (fp)
	{
		std::string contents;
		std::fseek(fp, 0, SEEK_END);
		contents.resize(std::ftell(fp));
		std::rewind(fp);

		// Stupid GCC bug. We do this to hide warnings.
		if(!std::fread(&contents[0], 1, contents.size(), fp))
			std::fclose(fp);
		else
			std::fclose(fp);

		return(contents);
	}
	throw(errno);
}

int main(int argc, char* argv[])
{
	mpi::environment env(argc, argv);
	mpi::communicator world;

	// Get requested number of walkers and make 
	// sure the number of processors is evenly divisible.
	int nwalkers = std::stoi(argv[1]);
	if(world.size() % nwalkers != 0)
	{
		if(world.rank() == 0)
			std::cerr << "The number of processors must be evenly "
			<< "divisible by the number of walkers." << std::endl;
		world.abort(-1);
	}

	// Determine walker ID using integer division and split 
	// communicator accordingly. 
	int wid = (int)world.rank() / nwalkers;
	auto walker = world.split(wid);

	// Silence of the lammps.
	char **largs = (char**) malloc(sizeof(char*) * 5);
	for(int i = 0; i < 5; ++i)
		largs[i] = (char*) malloc(sizeof(char) * 1024);
	sprintf(largs[0], " ");
	sprintf(largs[1], "-screen");
	sprintf(largs[2], "none");
	sprintf(largs[3], "-log");
	sprintf(largs[4], "log-MPI_ID-%d", wid);

	auto lammps = std::make_shared<LAMMPS>(5, largs, MPI_Comm(walker));

	// Read file on rank 0 proc, broadcast to all.
	std::string contents;
	if(world.rank() == 0)
		contents = GetFileContents(argv[2]);
	mpi::broadcast(world, contents, 0);

	// Go through lammps.
	std::string token;
	std::istringstream ss(contents);
	while(std::getline(ss, token, '\n'))
		lammps->input->one(token.c_str());

	// Initialize snapshot. 
	Snapshot snapshot(walker, wid);

	// Get hook from lammps modify.
	// Horrid, I know.
	auto fid = lammps->modify->find_fix("ssages");
	if(auto* hook = dynamic_cast<Hook*>(lammps->modify->fix[fid]))
	{
		hook->SetSnapshot(&snapshot);

		// Add methods and CV's here.
		///////Test Umbrella//////////////////////////////
		//hook->AddListener(new Umbrella({100.0}, {0}, 1));
		//hook->AddCV(new AtomCoordinateCV(1, 0));
		hook->AddListener(new MockMethod(1000));

		///////Test MetaDynamics//////////////////////////
		//hook->AddListener(new Meta(0.5, {0.05, 0.05}, 500, 1));
		//hook->AddCV(new AtomCoordinateCV(1, 0));
		//hook->AddCV(new AtomCoordinateCV(1, 1));
	}
	else
	{
		if(world.rank() == 0)
		{
			std::cerr << "Unable to dynamic cast hook. Error occurred" << std::endl;
			world.abort(-1);			
		}
	}

	// Run!
	std::string rline = "run " + std::string(argv[3]);
	lammps->input->one(rline.c_str());

	// Free.
	for(int i = 0; i < 5; ++i)
		free(largs[i]);
	free(largs);
}