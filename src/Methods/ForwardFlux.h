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

		enum Restart {NEW, LIBRARY, OLDCONFIG, NEWCONFIG, NONE};
		Restart _restart;

		// Output index file for storing information on where everything is.
		std::string _indexfilename; //User defined
		std::ofstream _indexfile;
		std::string _indexcontents;
		std::string _globalcontents;
		std::string _librarycontents;
		std::vector<Snapshot> _dumpconfigs;

		// Results file for end of simulation.
		std::string _resultsfilename; //User defined
		std::ofstream _resultsfile;
		std::string _resultscontents;

		// Location of the nodes to be used in determining FF interfaces
		std::vector<std::vector<double>> _centers;

		// Number of successes at a given interface
		std::vector<int> _successes;
		std::vector<int> _localsuccesses;
		std::vector<std::vector< int> > _paths;

		// Current interface FF is shooting from
		int _currentnode;

		// The current starting configuration that we are on
		int _currentstartingpoint;

		// User defined number of starting configs needed per walker before starting FF
		int _requiredconfigs;

		// Number that keeps track of configs this walker has generated
		int _currenthash;

		// Name of file of configuration where shooting from
		std::string _shootingconfigfile;

		// Number of shots each node takes per interface
		int _numshots;
		int _currentshot;

		// Flux
		int _fluxout;
		int _fluxin;

		// Probabilities
		std::vector<double> _weight;

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
				 int numshots,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _rd(), _gen(_rd()), _restart(LIBRARY),
		_indexfilename(indexfilename), _indexfile(), _indexcontents(""), _globalcontents(""), 
		_resultsfilename(resultsfilename), _resultsfile(), _resultscontents(""),
		_centers(centers), _successes(), _currentnode(currentinterface),
		_currentstartingpoint(-1), _requiredconfigs(requiredconfigs), _currenthash(10000*world.rank()), 
		_shootingconfigfile(""), _numshots(numshots), _currentshot(0), _fluxout(0), _fluxin(0)
		{
			_successes.resize(_centers.size());
			_localsuccesses.resize(_centers.size());
			for(size_t i = 0; i < _successes.size(); i++)
			{
				_successes[i] = 0;
				_localsuccesses[i] = 0;
			}
			_weight.resize(_requiredconfigs);

			if(newrun)
				_restart = NEW;
			else
				_restart = LIBRARY;
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
		// Return false if couldnt locate anything at a given interface
		bool ExtractInterfaceIndices(int interface, const std::string& contents,
									 std::vector<std::vector<std::string> >& InterfaceIndices);
		
		// Return the location of the nearest interface
		int AtInterface(const CVList& cvs)
		{
			std::vector<double> dists;
			dists.resize(_centers.size());

			// Record the difference between all cvs and all nodes
			for (size_t i = 0; i < _centers.size(); i++)
			{
				dists[i] = 0;
				for(size_t j = 0; j < cvs.size(); j++)
					dists[i]+=(cvs[j]->GetValue() - _centers[i][j])*(cvs[j]->GetValue() - _centers[i][j]);
			}

			return (std::min_element(dists.begin(), dists.end()) - dists.begin());
		}

		// store configuration (NEW in _libraryconfig, other in _dumpconfig)
		void StoreConfiguration(Snapshot* snapshot, int interface)
		{
			std::string dumpfilename = "dump_"+std::to_string(interface)+"_"+std::to_string(_currenthash)+".dump";
			_dumpconfigs.push_back(*snapshot);
			auto& snapshotID = _dumpconfigs.back().GetSnapshotID();
			auto& oldsnapshotID = snapshot->GetSnapshotID();
			
			snapshotID = dumpfilename;
			oldsnapshotID = dumpfilename;

			_currenthash++;
		}

		// Write out configuration file, this updates library and index file as well
		void WriteConfiguration(Snapshot* snapshot);

		// Read a given configuration from file and update snapshot
		void ReadConfiguration(Snapshot* snapshot, const std::string& dumpfilename);
		
		// Read a given configuration from _dumpcontents and update snapshot
		void ReadConfiguration(Snapshot* snapshot);

		// Clears everything to make way for a new run.
		void CleanUp();

		// Sets up new library for a new run because previous library had no full successes
		void SetUpNewLibrary(Snapshot* snapshot, const CVList& cvs);

		// Sets up the run from previous configuration
		bool SetUpRestartRun(Snapshot* snapshot);

		// Randomly picks a configuration from a given interface
		std::string PickConfiguration(int interface, const std::string& contents);

		// Randomly picks a configuration from list of _dumpconfigs
		void PickConfiguration();

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