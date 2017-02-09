/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Hythem Sidky <hsidky@nd.edu>
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
#include "ForwardFlux.h"
#include <iostream>
#include "../FileContents.h"
#include <random>

// This method involves a lot of bookkeeping. Typically the world node
// will hold gather all needed information and pass it along as it occurs.

namespace SSAGES
{
	void ForwardFlux::PreSimulation(Snapshot* /* snap */, const CVList& cvs)
	{
		if(cvs.size() > 1)
			throw BuildException({"Forwardflux currently only works with one cv."});

		switch(restart_)
		{
			case NEW:
			{
				indexfile_.open(indexfilename_.c_str(),std::ofstream::out | std::ofstream::trunc);
				resultsfile_.open(resultsfilename_.c_str(),std::ofstream::out | std::ofstream::trunc);
				libraryfile_.open(libraryfilename_.c_str(),std::ofstream::out | std::ofstream::trunc);
				break;
			}
			default:
			{
				indexfile_.open(indexfilename_.c_str(),std::ofstream::out | std::ofstream::app);
				resultsfile_.open(resultsfilename_.c_str(),std::ofstream::out | std::ofstream::app);
				libraryfile_.open(libraryfilename_.c_str(),std::ofstream::out | std::ofstream::app);
				break;
			}
		}

		currentnode_ = AtInterface(cvs);

		iteration_ = 0;
	}

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		switch(restart_)
		{
			case NEW:
			{
				SetUpNewLibrary(snapshot, cvs);
				return;
			}
			case LIBRARY:
			{
				if(world_.rank() == 0 && currentstartingpoint_ != 0)
				{
					for(auto& value : successes_)
						resultsfile_<<value<<" ";
					resultsfile_<<"\n";
				}

				mpi::all_reduce(world_, totalcontents_, globalcontents_, std::plus<std::string>());
				//Close local and global files
				if(world_.rank() == 0)
				{
					indexfile_<<globalcontents_<<std::endl;
				}
				totalcontents_ += indexcontents_;
				indexcontents_ = "";

				currentshot_ = 0;

				std::fill(localsuccesses_.begin(),localsuccesses_.end(),0);
				std::fill(successes_.begin(),successes_.end(),0);

				currentnode_ = 1;

				std::vector<std::vector<std::string> > tmp;
				if(!(ExtractInterfaceIndices(0, librarycontents_, tmp)))
				{
					std::cout<< "Could not locate any files";
					std::cout<< " at first interface!"<<std::endl;
					std::cout<< librarycontents_<<std::endl;
					PostSimulation(snapshot, cvs);
					world_.abort(0);
				}

				if(currentstartingpoint_ >= tmp.size())
				{
					if(world_.rank()==0)
						std::cout<<"Ending Simulation..."<<std::endl;
					PostSimulation(snapshot, cvs);
					world_.abort(0);
				}

				shootingconfigfile_ = tmp[currentstartingpoint_][1];
				if(world_.rank() == 0)
					currentconfig_ = GetFileContents(shootingconfigfile_.c_str());
				mpi::broadcast(world_, currentconfig_, 0);

				ReadConfiguration(snapshot, currentconfig_);
				restart_ = NONE;
				currentstartingpoint_++;
				return;
			}
			case NEWCONFIG:
			{
				currentshot_ = 0;
				mpi::all_reduce(world_, indexcontents_, globalcontents_, std::plus<std::string>());

				shootingconfigfile_ = PickConfiguration(currentnode_, globalcontents_);

				if(world_.rank() == 0)
					currentconfig_ = GetFileContents(shootingconfigfile_.c_str());
				
				mpi::broadcast(world_, currentconfig_, 0);
				ReadConfiguration(snapshot, currentconfig_);
				restart_ = NONE;

				return;
			}
			default:
			{
				break;
			}
		}

		// Locate the interface you are at and check if:
		// Returned to origin or at next interface
		unsigned int atinter = AtInterface(cvs);
		if(atinter == currentnode_ + 1 || atinter == 0)
		{
			currentnode_++;			
			// Check if you made it to the next one!
			if(atinter == currentnode_)
			{
				auto& ID = snapshot->GetSnapshotID();
				ID = "dump_"+std::to_string(currentnode_)+"_"+std::to_string(currenthash_)+".dump";
				currenthash_++;
				if(comm_.rank()==0)
				{
					WriteConfiguration(snapshot);
					localsuccesses_[currentnode_]++;
				}
			}

			if(currentshot_ < numshots_)
			{
				currentshot_++;
				ReadConfiguration(snapshot, currentconfig_);
				currentnode_--;
				iteration_++;
				return;
			}

			mpi::all_reduce(world_, localsuccesses_[currentnode_], successes_[currentnode_], std::plus<int>());
			if(successes_[currentnode_] > 0)
			{
				restart_ = NEWCONFIG;

				// If you have reached the final interface you are "done" with that path
				if(currentnode_ >= centers_.size()-1)
				{
					std::cout<< "Found finishing configuration! "<<std::endl;
					restart_ = LIBRARY;
				}
			}
			else
				restart_ = LIBRARY;

			iteration_++;
		}
	}

	void ForwardFlux::PostSimulation(Snapshot*, const CVList&)
	{
		//Close local and global files
		if(world_.rank() == 0)
		{
			resultsfile_<<"flux in: "<<fluxin_<<std::endl;
			resultsfile_<<"flux out: "<<fluxout_<<std::endl;

			indexfile_.close();
			resultsfile_.close();
			libraryfile_.close();
		}
	}

	// Setting up new run, so setup new starting configurations at the first interface
	void ForwardFlux::SetUpNewLibrary(Snapshot* snapshot, const CVList& cvs)
	{
		// Get CV values check if at next interface, if so store configuration
		int interface = AtInterface(cvs);
		shootingconfigfile_ = "Origin";

		// Flux out of A
		if(interface == 1 && currentnode_ == 0)
		{
			auto& ID = snapshot->GetSnapshotID();
			ID = "dump_"+std::to_string(interface)+"_"+std::to_string(currenthash_)+".dump";
			currenthash_++;
			if(comm_.rank()==0)
			{
				WriteConfiguration(snapshot);
				localsuccesses_[currentnode_]++;
			}

			fluxout_++;
		}
		// Flux back in towards A
		else if(interface == 0 && currentnode_ == 1)
			fluxin_++;

		currentnode_ = interface;

		std::vector<std::vector<std::string> > TempLibrary;
		ExtractInterfaceIndices(0, indexcontents_, TempLibrary);

		if(TempLibrary.size() >= requiredconfigs_ && currentnode_ == 0)
		{
			restart_ = LIBRARY;
			mpi::all_reduce(world_, indexcontents_, librarycontents_, std::plus<std::string>());
			currentstartingpoint_ = 0;
			if(world_.rank()==0)
				libraryfile_<<librarycontents_<<std::endl;
		}
	}

	void ForwardFlux::WriteConfiguration(Snapshot* snapshot)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();
		const auto& dumpfilename = snapshot->GetSnapshotID();

		// Write the dump file out
		std::ofstream dumpfile;
 		dumpfile.open(dumpfilename.c_str());

 		for(size_t i = 0; i< atomID.size(); i++)
 		{
 			dumpfile<<atomID[i]<<" ";
 			dumpfile<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" ";
 			dumpfile<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
		}

		std::vector<std::string> tmpstr;
		tmpstr.push_back(std::to_string(currentnode_));
		tmpstr.push_back(dumpfilename);
 		tmpstr.push_back(shootingconfigfile_);

 		// Update index file of new configuration
 		indexcontents_ += tmpstr[0]+" "+tmpstr[1]+" "+tmpstr[2]+"\n";
 		dumpfile.close();
	}

	// Extract all indices for a given interface and contetns. 
	// Return false if couldnt locate anything at a given interface
	bool ForwardFlux::ExtractInterfaceIndices(unsigned int interface, const std::string& contents,
											 std::vector<std::vector<std::string> >& InterfaceIndices)
	{
		//Extract configuration indices for a given interface int
		std::istringstream f(contents);
		std::string line;
		while (std::getline(f, line))
		{
			std::string buf; // Have a buffer string
			std::stringstream ss(line); // Insert the string into a stream
			std::vector<std::string> tokens; // Create vector to hold our words

			while (ss >> buf)
			    tokens.push_back(buf);

			if(std::stoul(tokens[0]) == interface)
				InterfaceIndices.push_back(tokens);
		}

		if(InterfaceIndices.size() == 0)
			return false;

		return true;
	}

	void ForwardFlux::ReadConfiguration(Snapshot* snapshot, const std::string& FileContents)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		auto& ID = snapshot->GetSnapshotID();

		ID = shootingconfigfile_;

		//Extract currentconfig information
		std::istringstream f(FileContents);
		std::string line;
		while (std::getline(f, line))
		{
			int atomindex = -1;
			std::string buf; // Have a buffer string
			std::stringstream ss(line); // Insert the string into a stream
			std::vector<std::string> tokens; // Create vector to hold our words

			while (ss >> buf)
			    tokens.push_back(buf);

			if(tokens.size() != 7)
			{
				std::cout<<"error, incorrect line format in "<<shootingconfigfile_<<" on line: "<<std::endl;
				std::cout<<line<<std::endl;
				world_.abort(-1);	
			}

			for(size_t i=0; i < atomID.size(); i++)
			{
				if(atomID[i] == std::stoi(tokens[0]))
					atomindex = i;
			}

			if(atomindex < 0)
			{
				std::cout<<"error, could not locate atomID "<<tokens[0]<<" from dumpfile"<<std::endl;
				world_.abort(-1);
			}

			positions[atomindex][0] = std::stod(tokens[1]);
			positions[atomindex][1] = std::stod(tokens[2]);
			positions[atomindex][2] = std::stod(tokens[3]);
			velocities[atomindex][0] = std::stod(tokens[4]);
			velocities[atomindex][1] = std::stod(tokens[5]);
			velocities[atomindex][2] = std::stod(tokens[6]);

			for(auto& force : forces)
				force.setZero();
		}
	}

	// Pick a random configuration
	std::string ForwardFlux::PickConfiguration(unsigned int interface, const std::string& contents)
	{
		std::string configfilename;
		if(world_.rank() == 0)
		{
			std::vector<std::vector<std::string> > files; 
			if(!(ExtractInterfaceIndices(interface, contents, files)))
			{
				std::cout<< "Could not locate any files at interface ";
				std::cout<< interface<<" in PickCOnfiguration!"<<std::endl;
				std::cout<< contents<<std::endl;
				world_.abort(-1);
			}

			std::uniform_int_distribution<> dis(0, files.size()-1);
			int configfile = dis(gen_);
			configfilename = files[configfile][1];
		}

		mpi::broadcast(world_, configfilename, 0);

		return configfilename;
	}
}

