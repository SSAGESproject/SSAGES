/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Joshua Lequieu <lequieu@uchicago.edu>
 *                Hadi Ramezani-Dakhel <ramezani@uchicago.edu>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <random>
#include "../FileContents.h"
#include "../Drivers/DriverException.h"

namespace mpi = boost::mpi;
namespace SSAGES
{
	//! ForwardFlux sampling method
	/*!
	 * \ingroup Methods
     * The notation used here is drawn largely from Allen, Valeriani and Rein ten Wolde. J. Phys.: Condens. Matter (2009) 21:463102. 
     * We recommend referring to this review if the reader is unfamiliar with the method, or our variable naming conventions.
	 */
	class ForwardFlux : public Method
	{
	private:

		//! Number of FFS interfaces
		std::vector<double> _ninterfaces;

        //! Current Interface
        unsigned int _currentinterface;

		//! Previous cv position, used to determine if you've crossed an interface since last time
        double _cvposition_previous;

		//! Number of configurations to collect at lambda0 (first interface) 
		unsigned int _N0 ;

        //! Data structure that holds a Library N0 configurations at lambda0
        std::vector<FFSConfiguration> Lambda0ConfigLibrary

        //! Total Simulation Time spent in accumulating \ _N0
        double _N0TotalSimTime;

        //! Flux of trajectories out of state A. Denoted PhiA0 over h_A in Allen2009.
        double _fluxA0;

        //! Number of trials to attemts from each interface
        std::vector<double> _M;

        //! Flag to determine wheter fluxA0 should be calculated
        bool _computefluxA0;

        //! Stores what 'mode' of FFS we're in. 
        /*!
         *  Options:
         *   - computefluxA0
         *   - WHAT OTHERS?
         */
        unsigned int _FFSmode;

        struct FFSConfigID
        {
            unsigned int lambda; //!< Interface number
            unsigned int n;      //!< Configuration Number
            unsigned int a;      //!< Attempt number
            FFSConfigID previous; //!< ID of FFSConfiguration that I came from
        };

        //! The current FFSConfigID of this MPI process
        FFSConfiguration myFFSConfigID;

        //! Data structure that holds FFSConfigurations
        /*!
         *  When a given processor reaches an interface, it pulls a config from this Queue to figure out what it should do next
         *  This object should be syncronized between all FFS walkers (is walker the correct terminology here?)
         */
        std::queue<FFSConfigID> FFSConfigIDQueue; 
 
        
        //! Function checks if configuration has returned to A
        bool HasReturnedToA();

        //! Function checks if configuration has crossed interface specified since the last check
        bool HasCrossedInterface(unsigned int interface);

        //! Function checks if FFS is Finished, returns bool with result
        //! See if interface is the last one, and the queue is empty, etc
        bool CheckIfFinishedMethod();

        //! Write a file corresponding to FFSConfigID from current snapshot
        bool WriteFFSConfiguration(Snapshot *,FFSConfigID);

        //! Read a file corresponding to a FFSConfigID into current snapshot
        bool ReadFFSConfiguration(Snapshot *,FFSConfigID);

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param indexfilename File name of index file.
		 * \param libraryfilename File name for the library file.
		 * \param resultsfilename File name of results file.
		 * \param centers List of centers.
		 * \param requiredconfigs Number of required configurations.
		 * \param numshots Number of shots each node takes.
		 * \param frequency Frequency with which this method is invoked.
		 *
		 * Create instance of Forward Flux with centers "centers".
		 */
		ForwardFlux(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 std::string indexfilename,
				 std::string libraryfilename,
				 std::string resultsfilename,
				 std::vector<double> centers,
				 unsigned int requiredconfigs,
				 int numshots,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _rd(), _gen(_rd()), _restart(NEW),
		_indexfilename(indexfilename), _indexcontents(""), _globalcontents(""),
		_totalcontents(""), _libraryfilename(libraryfilename), _librarycontents(""),
		_resultsfilename(resultsfilename), _resultscontents(""), 
		_centers(centers), _currentnode(0), _currentstartingpoint(0), 
		_requiredconfigs(requiredconfigs), _currenthash(1000000*world.rank()), 
		_shootingconfigfile(""), _numshots(numshots), _currentshot(0),
		_fluxout(0), _fluxin(0)
		{
			_successes.resize(_centers.size());
			_localsuccesses.resize(_centers.size());
			for(size_t i = 0; i < _successes.size(); i++)
			{
				_successes[i] = 0;
				_localsuccesses[i] = 0;
			}
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

		//! Extract all indices for a given interface.
		/*!
		 * \param interface Index of the interface.
		 * \param contents The contents will be written to this string.
		 * \param InterfaceIndices List of interface indices.
		 * \return \c False if nothing could be located at a given interface.
		 */
		bool ExtractInterfaceIndices(unsigned int interface, const std::string& contents,
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
				dists[i] = (cvs[0]->GetValue() - _centers[i])*(cvs[0]->GetValue() - _centers[i]);

			return (std::min_element(dists.begin(), dists.end()) - dists.begin());
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

		//! Sets up new library for a new run because previous library had no full successes
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void SetUpNewLibrary(Snapshot* snapshot, const CVList& cvs);

		//! Randomly picks a configuration from a given interface.
		/*!
		 * \param interface Index of the interface.
		 * \param contents String containing the configuration.
		 * \return Configuration file name.
		 */
		std::string PickConfiguration(unsigned int interface, const std::string& contents);

		//! Set restart type.
		/*!
		 * \param restart String specifying the restart type.
		 *
		 * Set the type of restart. Possible values are:
		 *
		 * String Value     | Behavior
		 * ------------     | --------
		 * "new_library"    | Start method from scratch.
		 * "from_library"   | Load a previous library and continue run.
		 * "from_interface" | Start with a new configuration.
		 * "none"           | No restart.
		 *
		 * If a value different from one of the above is given, the method will
		 * throw a BuildException.
		 */
		void SetRestart(std::string restart)
		{
			if(restart == "new_library")
				_restart = NEW;
			else if (restart == "from_library")
				_restart = LIBRARY;
			else if (restart == "from_interface")
				_restart = NEWCONFIG;
			else if (restart == "none")
				_restart = NONE;
			else
				throw BuildException({"Unknown restart type for Forward Flux Method!"});
		}

		//! Set starting point of the library.
		/*!
		 * \param lpoint New starting point for the library.
		 */
		void SetLibraryPoint(unsigned int lpoint){_currentstartingpoint = lpoint;}

		//! Set hash value.
		/*!
		 * \param hash New hash value.
		 */
		void SetHash(int hash){_currenthash = hash;}

		//! Set contents for the index.
		/*!
		 * \param contents New string for the index contents.
		 */
		void SetIndexContents(std::string contents){_indexcontents = contents;}

		//! Set the value for the current shot.
		/*!
		 * \param atshot New value for the current shot.
		 */
		void SetAtShot(int atshot){_currentshot = atshot;}

		//! Set the number of successes for the interfaces.
		/*!
		 * \param succ Vector storing the number of successes for all interfaces.
		 *
		 * Note that the size of the vector \c succ must equal the number of
		 * interfaces in the method.
		 */
		void SetSuccesses(std::vector<int> succ)
		{
			if(_localsuccesses.size() != succ.size())
				throw BuildException({"Number of interfaces does not match local successes."});

			for(size_t i = 0; i<succ.size(); i++)
				_localsuccesses[i] = succ[i];
		}

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override
		{
			//Needed to run
			json["type"] = "ForwardFlux";
			json["index_file"] = _indexfilename;
			json["library_file"] = _libraryfilename;
			json["results_file"] = _resultsfilename;
			json["generate_configs"] = _requiredconfigs;
			json["shots"] = _numshots;
			for(auto& c : _centers)
				json["centers"].append(c);
			
			// Needed for Restart:
			std::string restart;
			switch(_restart)
			{
				case NEW: 
				{
					restart = "new_library";
					break;
				}
				case LIBRARY: 
				{
					restart = "from_library";
					break;
				}
				case NEWCONFIG: 
				{
					restart = "from_interface";
					break;
				}
				default: 
				{
					restart = "none";
				}
			}

			json["restart_type"] = restart;
			json["library_point"] = _currentstartingpoint;
			json["current_hash"] = _currenthash;
			json["index_contents"] = _indexcontents;
			for(auto s : _localsuccesses)
				json["successes"].append(s);

			json["current_shot"] = _currentshot;

		}

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
