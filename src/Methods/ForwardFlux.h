#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include "../FileContents.h"

namespace mpi = boost::mpi;
namespace SSAGES
{
	// ForwardFlux sampling method 
	class ForwardFlux : public Method
	{
	private:

		std::random_device _rd;
		std::mt19937 _gen;

		// Output index file for storing information on where everything is.
		std::string _indexfilename; //User defined
		std::ofstream _indexfile;
		std::string _indexcontents;

		// Results file for end of simulation.
		std::string _resultsfilename; //User defined
		std::ofstream _resultsfile;
		std::string _resultscontents;

		// Location of the nodes to be used in determining FF interfaces
		std::vector<std::vector<double>> _centers;

		// Number of successes at a given interface
		std::vector<int> _successes;

		// User defined if we need to create a library of new configs or not
		bool _newrun;
		bool _restartfromlibrary;
		bool _restartfrominterface;

		// Current interface FF is shooting from
		int _currentinterface;

		// The current starting configuration that we are on
		int _currentstartingpoint;

		// User defined number of starting configs needed per walker before starting FF
		int _requiredconfigs;

		// Number that keeps track of configs this walker has generated
		int _currenthash;

		// Name of file of configuration where shooting from
		std::string _shootingconfigfile;

		// Flux
		int _fluxout;
		int _fluxin;

	public:
		// Create instance of Forward Flux with centers "centers". 
		ForwardFlux(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 std::string indexfilename,
				 std::string resultsfilename,
				 int currentinterface,
				 std::vector<std::vector<double> > centers,
				 bool newrun,
				 int requiredconfigs,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _rd(), _gen(_rd()), _indexfilename(indexfilename),
		_indexfile(), _indexcontents(), _resultsfilename(resultsfilename), _resultsfile(),
		_resultscontents(), _centers(centers), _successes(), _newrun(newrun), _restartfromlibrary(),
		_restartfrominterface(), _currentinterface(currentinterface),_currentstartingpoint(),
		_requiredconfigs(requiredconfigs), _currenthash(), _shootingconfigfile(),_fluxout(0), _fluxin(0)
		{
			_successes.resize(_centers.size());
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;
		
		void Serialize(Json::Value& json) const override
		{

		}

		// Extract all indices for a given interface. 
		// Return true if couldnt locate anything at a given interface
		bool ExtractInterfaceIndices(int interface, std::vector<std::vector<std::string> >& InterfaceIndices);
		
		// Return the location of the nearest interface
		int AtInterface(const CVList& cvs);

		// Write out configuration file, this updates library and index file as well
		void WriteConfiguration(Snapshot* snapshot);

		// Read a given configuration and update snapshot
		void ReadConfiguration(Snapshot* snapshot, std::string dumpfilename);

		// Clears everything to make way for a new run. Will overwrite and clear files as well
		void ClearFiles();

		// Sets up new library for a new run because previous library had no full successes
		void SetUpNewRun(Snapshot* snapshot, const CVList& cvs);

		// Randomly picks a configuration from a given interface
		std::string PickConfiguration(int interface);

	};
}


/*
File Formats:
_indexfile
interface(some integer) dump_file_name(a string that contains interface and trial number)
example: 1 dump_1_10.xyz

dumpfile
atomid posx posy posz vx vy vz


*/