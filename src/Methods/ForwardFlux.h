/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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
	 */
	class ForwardFlux : public Method
	{
	private:

		std::random_device rd_; //!< Random number generator.
		std::mt19937 gen_; //!< Alternative random number generator.

		//! Possible restart values.
		enum Restart {NEW, LIBRARY, NEWCONFIG, NONE};

		//! Restart value.
		Restart restart_;

		// Output index file for storing information on where everything is.
		std::string indexfilename_; //!< File name for index file.
		std::ofstream indexfile_; //!< File stream for index file.
		std::string indexcontents_; //!< Contents of index file.
		std::string globalcontents_; //!< Global contents.
		std::string totalcontents_; //!< On going record of all paths.

		std::string libraryfilename_; //!< File name for index file.
		std::ofstream libraryfile_; //!< File stream for index file.
		std::string librarycontents_; //!< Library contents.
		
		std::string currentconfig_; //!< Current configuration.

		// Results file for end of simulation.
		std::string resultsfilename_; //!< File name for simulation results.
		std::ofstream resultsfile_; //!< File stream for simulation results.
		std::string resultscontents_; //!< Content of simulation results.

		//! Location of the nodes to be used in determining FF interfaces.
		std::vector<double> centers_;

		//! Number of successes at a given interface.
		std::vector<int> successes_;

		//! Number of local successes.
		std::vector<int> localsuccesses_;

		//! Current interface FF is shooting from.
		unsigned int currentnode_;

		//! The current starting configuration that we are on.
		unsigned int currentstartingpoint_;

		//! User defined number of starting configs needed per walker before starting FF.
		unsigned int requiredconfigs_;

		//! Number that keeps track of configs this walker has generated.
		int currenthash_;

		//! Name of file of configuration where shooting from.
		std::string shootingconfigfile_;

		//! Number of shots each node takes per interface
		int numshots_;

		//! Index of the current shot.
		int currentshot_;

		// Flux
		int fluxout_; //!< Flux out.
		int fluxin_; //!< Flux in.

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
		Method(frequency, world, comm), rd_(), gen_(rd_()), restart_(NEW),
		indexfilename_(indexfilename), indexcontents_(""), globalcontents_(""),
		totalcontents_(""), libraryfilename_(libraryfilename), librarycontents_(""),
		resultsfilename_(resultsfilename), resultscontents_(""), 
		centers_(centers), currentnode_(0), currentstartingpoint_(0), 
		requiredconfigs_(requiredconfigs), currenthash_(1000000*world.rank()), 
		shootingconfigfile_(""), numshots_(numshots), currentshot_(0),
		fluxout_(0), fluxin_(0)
		{
			successes_.resize(centers_.size());
			localsuccesses_.resize(centers_.size());
			for(size_t i = 0; i < successes_.size(); i++)
			{
				successes_[i] = 0;
				localsuccesses_[i] = 0;
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
			dists.resize(centers_.size());

			// Record the difference between all cvs and all nodes
			for (size_t i = 0; i < centers_.size(); i++)
				dists[i] = (cvs[0]->GetValue() - centers_[i])*(cvs[0]->GetValue() - centers_[i]);

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
				restart_ = NEW;
			else if (restart == "from_library")
				restart_ = LIBRARY;
			else if (restart == "from_interface")
				restart_ = NEWCONFIG;
			else if (restart == "none")
				restart_ = NONE;
			else
				throw BuildException({"Unknown restart type for Forward Flux Method!"});
		}

		//! Set starting point of the library.
		/*!
		 * \param lpoint New starting point for the library.
		 */
		void SetLibraryPoint(unsigned int lpoint){currentstartingpoint_ = lpoint;}

		//! Set hash value.
		/*!
		 * \param hash New hash value.
		 */
		void SetHash(int hash){currenthash_ = hash;}

		//! Set contents for the index.
		/*!
		 * \param contents New string for the index contents.
		 */
		void SetIndexContents(std::string contents){indexcontents_ = contents;}

		//! Set the value for the current shot.
		/*!
		 * \param atshot New value for the current shot.
		 */
		void SetAtShot(int atshot){currentshot_ = atshot;}

		//! Set the number of successes for the interfaces.
		/*!
		 * \param succ Vector storing the number of successes for all interfaces.
		 *
		 * Note that the size of the vector \c succ must equal the number of
		 * interfaces in the method.
		 */
		void SetSuccesses(std::vector<int> succ)
		{
			if(localsuccesses_.size() != succ.size())
				throw BuildException({"Number of interfaces does not match local successes."});

			for(size_t i = 0; i<succ.size(); i++)
				localsuccesses_[i] = succ[i];
		}

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override
		{
			//Needed to run
			json["type"] = "ForwardFlux";
			json["index_file"] = indexfilename_;
			json["library_file"] = libraryfilename_;
			json["results_file"] = resultsfilename_;
			json["generate_configs"] = requiredconfigs_;
			json["shots"] = numshots_;
			for(auto& c : centers_)
				json["centers"].append(c);
			
			// Needed for Restart:
			std::string restart;
			switch(restart_)
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
			json["library_point"] = currentstartingpoint_;
			json["current_hash"] = currenthash_;
			json["index_contents"] = indexcontents_;
			for(auto s : localsuccesses_)
				json["successes"].append(s);

			json["current_shot"] = currentshot_;

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
