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


    }

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
        //check if we want to check FFS interfaces this timestep

        // if _computefluxA0
        if (_FFSmode == "InitialFlux" ){
            ComputeInitialFlux(); 
        }
        // Else normal forward flux
        else if (_FFSmode == "ForwardFlux" ){
          //check if I've crossed the next interface
          
          //if so update some relevant quantities and mpi them across procs

          //assign myFFSConfigID = NULL

          // reassign configs that reached interface a new FFSConfigID
           
        }
        // Other modes?

    }

	void ForwardFlux::PostSimulation(Snapshot*, const CVList&)
	{
		
	}

    bool ForwardFlux::HasReturnedToA(Snapshot* snapshot){


    }

	void ForwardFlux::ComputeInitialFlux(){
        //check if I've crossed the first interface (lambda 0)

        //need to sync variables between processors
    }

	void ForwardFlux::ComputeTransitionProbabilities(){

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
		tmpstr.push_back(std::to_string(_currentnode));
		tmpstr.push_back(dumpfilename);
 		tmpstr.push_back(_shootingconfigfile);

 		// Update index file of new configuration
 		_indexcontents += tmpstr[0]+" "+tmpstr[1]+" "+tmpstr[2]+"\n";
 		dumpfile.close();
	}


	void ForwardFlux::ReadConfiguration(Snapshot* snapshot, const std::string& FileContents)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		auto& ID = snapshot->GetSnapshotID();

		ID = _shootingconfigfile;

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
				std::cout<<"error, incorrect line format in "<<_shootingconfigfile<<" on line: "<<std::endl;
				std::cout<<line<<std::endl;
				_world.abort(-1);	
			}

			for(size_t i=0; i < atomID.size(); i++)
			{
				if(atomID[i] == std::stoi(tokens[0]))
					atomindex = i;
			}

			if(atomindex < 0)
			{
				std::cout<<"error, could not locate atomID "<<tokens[0]<<" from dumpfile"<<std::endl;
				_world.abort(-1);
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
		if(_world.rank() == 0)
		{
			std::vector<std::vector<std::string> > files; 
			if(!(ExtractInterfaceIndices(interface, contents, files)))
			{
				std::cout<< "Could not locate any files at interface ";
				std::cout<< interface<<" in PickCOnfiguration!"<<std::endl;
				std::cout<< contents<<std::endl;
				_world.abort(-1);
			}

			std::uniform_int_distribution<> dis(0, files.size()-1);
			int configfile = dis(_gen);
			configfilename = files[configfile][1];
		}

		mpi::broadcast(_world, configfilename, 0);

		return configfilename;
	}
}

