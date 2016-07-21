#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <random>
#include "../FileContents.h"

namespace mpi = boost::mpi;
namespace SSAGES
{
	//! ForwardFlux sampling method
	/*!
	 * \ingroup Methods
	 */
	class ForwardFlux : public Method
	{
	private:

		std::random_device _rd; //!< Random number generator.
		std::mt19937 _gen; //!< Alternative random number generator.

		//! Possible restart values.
		enum Restart {NEW, LIBRARY, OLDCONFIG, NEWCONFIG, RESTART, NONE};

		//! Restart value.
		Restart _restart;

		// Output index file for storing information on where everything is.
		std::string _indexfilename; //!< File name for index file.
		std::ofstream _indexfile; //!< File stream for index file.
		std::string _indexcontents; //!< Contents of index file.
		std::string _globalcontents; //!< Global contents.
		std::string _librarycontents; //!< Library contents.
		std::string _restartfilename; //!< File name for restart file.
		std::string _currentconfig; //!< Current configuration.
		std::vector<Snapshot> _dumpconfigs; //!< Configurations to write out.

		// Results file for end of simulation.
		std::string _resultsfilename; //!< File name for simulation results.
		std::ofstream _resultsfile; //!< File stream for simulation results.
		std::string _resultscontents; //!< Content of simulation results.

		//! Location of the nodes to be used in determining FF interfaces.
		std::vector<std::vector<double>> _centers;

		//! Number of successes at a given interface.
		std::vector<int> _successes;

		//! Number of local successes.
		std::vector<int> _localsuccesses;

		//! List of paths.
		std::vector<std::vector< int> > _paths;

		//! Current interface FF is shooting from.
		int _currentnode;

		//! The current starting configuration that we are on.
		int _currentstartingpoint;

		//! User defined number of starting configs needed per walker before starting FF.
		int _requiredconfigs;

		//! Number that keeps track of configs this walker has generated.
		int _currenthash;

		//! Name of file of configuration where shooting from.
		std::string _shootingconfigfile;

		//! Number of shots each node takes per interface
		int _numshots;

		//! Index of the current shot.
		int _currentshot;

		// Flux
		int _fluxout; //!< Flux out.
		int _fluxin; //!< Flux in.

		//! Probabilities
		std::vector<double> _weight;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param indexfilename File name of index file.
		 * \param resultsfilename File name of results file.
		 * \param restartfilename File name for restart file.
		 * \param currentinterface Index of the current interface.
		 * \param centers List of centers.
		 * \param newrun If \c True start a new run.
		 * \param requiredconfigs Number of required configurations.
		 * \param numshots Number of shots each node takes.
		 * \param frequency Frequency with which this method is invoked.
		 *
		 * Create instance of Forward Flux with centers "centers".
		 */
		ForwardFlux(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 std::string indexfilename,
				 std::string resultsfilename,
				 std::string restartfilename,
				 int currentinterface,
				 std::vector<std::vector<double> > centers,
				 bool newrun,
				 int requiredconfigs,
				 int numshots,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _rd(), _gen(_rd()), _restart(LIBRARY),
		_indexfilename(indexfilename), _indexfile(), _indexcontents(""), _globalcontents(""), 
		_resultsfilename(resultsfilename), _restartfilename(restartfilename),
		_resultsfile(), _resultscontents(""), _centers(centers), _successes(), _currentnode(currentinterface),
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
			else if(_restartfilename != "none")
				_restart = RESTART;
			else
				_restart = LIBRARY;
		}

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;
		
		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning Serialization not implemented yet!
		 */
		void Serialize(Json::Value& json) const override
		{

		}

		//! Extract all indices for a given interface.
		/*!
		 * \param interface Index of the interface.
		 * \param contents The contents will be written to this string.
		 * \param InterfaceIndices List of interface indices.
		 * \return \c False if nothing could be located at a given interface.
		 */
		bool ExtractInterfaceIndices(int interface, const std::string& contents,
									 std::vector<std::vector<std::string> >& InterfaceIndices);
		
		//! Return the location of the nearest interface
		/*!
		 * \param cvs List of CVs.
		 * \return Location of nearest interface.
		 */
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

		//! Store configuration
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param interface Index of the interface.
		 *
		 * NEW: store in _libraryconfig, other in _dumpconfig
		 */
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

		//! Write out configuration file
		/*!
		 * \param snapshot Current simulation snapshot.
		 *
		 * This updates library and index file as well.
		 */
		void WriteConfiguration(Snapshot* snapshot);

		//! Read a given configuration from file and update snapshot
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param dumpfilename File name of the dump file.
		 */
		void ReadConfiguration(Snapshot* snapshot, const std::string& dumpfilename);
		
		//! Read a given configuration from _dumpcontents and update snapshot
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void ReadConfiguration(Snapshot* snapshot);

		//! Clears everything to make way for a new run.
		void CleanUp();

		//! Sets up new library for a new run because previous library had no full successes
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void SetUpNewLibrary(Snapshot* snapshot, const CVList& cvs);

		//! Sets up the run from previous configuration
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \return \c False if the run could not be set up successfully.
		 */
		bool SetUpRestartRun(Snapshot* snapshot);

		//! Randomly picks a configuration from a given interface
		/*!
		 * \param interface Index of the interface.
		 * \param contents String containing the configuration.
		 * \return Configuration file name.
		 */
		std::string PickConfiguration(int interface, const std::string& contents);

		//! Randomly picks a configuration from list of _dumpconfigs
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
