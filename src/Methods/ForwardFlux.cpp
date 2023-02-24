/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Hythem Sidky <hsidky@nd.edu>
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
#include "ForwardFlux.h"
#include "DirectForwardFlux.h"
#include "FileContents.h"
#include "CVs/CVManager.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"
#include "sys/stat.h"
#include <iostream>
#include <queue>

using namespace Json;

namespace SSAGES
{
	void ForwardFlux::PreSimulation(Snapshot* /* snap */, const CVManager& cvmanager)
	{
		auto cvs = cvmanager.GetCVs(cvmask_);
		if(cvs.size() > 1)
			throw BuildException({"Forwardflux currently only works with one cv."});

        std::cout << "\nWARNING! MAKE SURE LAMMPS GIVES A DIFFERENT RANDOM SEED TO EACH PROCESSOR, OTHERWISE EACH FFS TRAJ WILL BE IDENTICAL!\n";        
    }

	void ForwardFlux::PostSimulation(Snapshot* snapshot, const CVManager&)
	{
        std::cout << "Post simulation\n";

        ComputeCommittorProbability(snapshot);
        ComputeTransitionProbabilities();

        if (_saveTrajectories)
        	ReconstructTrajectories(snapshot);	
        
		MPI_Abort(world_, EXIT_FAILURE); //more elegant solution?
	}

	void ForwardFlux::CheckInitialStructure(const CVList& cvs)
	{
		if (_initialFluxFlag)
		{
			std::cout << "Running initial Flux calculations" << std::endl;

			// Check if we are in State A
			_cvvalue = cvs[0]->GetValue();
			_cvvalue_previous = _cvvalue;
			double _firstInterfaceLocation = _interfaces[0];
			if ( _cvvalue > _firstInterfaceLocation) 
			{
				std::cerr << "Please provide an initial configuration in State A. Exiting ...." << std::endl;
				MPI_Abort(world_, EXIT_FAILURE);
			}
		}
	}
    
    int ForwardFlux::HasCrossedInterface(double current, double prev, unsigned int i)
	{
    	double interface_location = _interfaces[i];
        if (_interfaces_increase)
		{
          if ((prev <= interface_location) && (current >= interface_location))
              return 1;
          else if ((prev >= interface_location) && (current <= interface_location))
              return -1;
          else
              return 0;
        }
        else
		{
          if ((prev >= interface_location) && (current <= interface_location))
              return 1;
          else if ((prev <= interface_location) && (current >= interface_location))
              return -1;
          else
              return 0;
        }
    }

    bool ForwardFlux::HasReturnedToA(double current)
	{
        double interface_location = _interfaces[0];
        if (_interfaces_increase)
		{
        	if (current < interface_location) return true;
        	else return false;
        }
        else
		{
        	if (current > interface_location) return true;
        	else return false;
        }
    }

	void ForwardFlux::ComputeInitialFlux(Snapshot* snapshot, const CVList& cvs)
	{
        _cvvalue = cvs[0]->GetValue();

        //check if we've crossed the first interface (lambda 0)
        int hascrossed = HasCrossedInterface(_cvvalue, _cvvalue_previous, 0);
        unsigned int success_local = false;
        std::vector<unsigned int> successes (world_.size(),0);

        if (hascrossed == 1)
        	success_local = true;

        //for each traj that crossed to lambda0 in forward direction, we need to write it to disk (FFSConfigurationFile)
        MPI_Allgather(&success_local,1,MPI_UNSIGNED,successes.data(),1,MPI_UNSIGNED,world_);

        int success_count = 0;
        for (int i = 0; i < world_.size(); i++)
		{
			// Since we are in State A, the values of lprev, nprev, aprev are all zero.
			if (successes[i] == true)
			{ 
				if (i == world_.rank())
				{
					int l,n,a,lprev,nprev,aprev;
					//update ffsconfigid's l,n,a
					//note lprev == l for initial interface
					l = 0;
					n = _N[0] + success_count;
					a = 0;
					lprev = l;
					nprev = n;
					aprev = a;

					FFSConfigID newid = FFSConfigID(l,n,a,lprev,nprev,aprev);
					Lambda0ConfigLibrary.emplace_back(l,n,a,lprev,nprev,aprev);
					WriteFFSConfiguration(snapshot,newid,1);
				}
				success_count++;
			}    
        }  

        // all procs update correctly
        _N[0] += success_count;

        // If not in B, increment the time
        double N0SimTime_local = 0;
        double N0SimTime;
        int reachedB = HasCrossedInterface(_cvvalue, _cvvalue_previous, _ninterfaces-1);
        // Question: This condition says that if we cross the last interface-> stop counting. We need to stop counting as long as we are in state B.
        if (!(reachedB == 1) && (_cvvalue < _interfaces.back()))
        	N0SimTime_local++;  //or += frequency if FFS isn't called on every step!

        // Allreduce then increment total
        MPI_Allreduce(&N0SimTime_local, &N0SimTime, 1, MPI_DOUBLE, MPI_SUM,world_);
        _N0TotalSimTime += N0SimTime;

        //print some info
        if (success_local)
		{
        	std::cout << "Iteration: "<< iteration_ << ", proc " << world_.rank() << std::endl;
            std::cout << "Successful attempt. (cvvalue_previous: " << _cvvalue_previous << " cvvalue " << _cvvalue << " )" << std::endl;
            std::cout << "# of successes:               " << _N[0] << std::endl;
            std::cout << "required # of configurations: " << _N0Target << std::endl;
        }

        // Check if the required number of initial configurations are created, if so print a message, and compute the initial flux.
        if (_N[0] >= _N0Target)
		{
            std::cout << "Initial flux calculation was successfully completed" << std::endl;
            // Call a function to compute the actual flux and output the data
            WriteInitialFlux();
            //responsible for setting _initialFluxFlag = false when finished
            _initialFluxFlag = false;
        } 

        _cvvalue_previous = _cvvalue;
        iteration_++;
    }

	void ForwardFlux::WriteInitialFlux()
	{
		std::ofstream file;
		std::string filename = _output_directory + "/initial_flux_value.dat";
		file.open(filename.c_str());
		if(!file)
		{
			std::cerr << "Error! Unable to write " << filename << std::endl;
			MPI_Abort(world_, EXIT_FAILURE);
		}
		_fluxA0 = (double) (_N[0] / _N0TotalSimTime);
		file << "number of walkers: " << world_.size() << std::endl;
		file << "number of iterations: " << iteration_ << std::endl;
		file << "Total simulation time: " << _N0TotalSimTime << std::endl;
		file << "Initial flux: " << _fluxA0 << std::endl;
		file.close();
	}

	void ForwardFlux::ComputeTransitionProbabilities()
	{
        double Ptotal = 1;
        for (unsigned int i = 0; i < _ninterfaces -1; i++)
            Ptotal *= _P[i];

        _rate = Ptotal*_fluxA0;
        
        //write file
        std::ofstream file;
        std::string filename = _output_directory + "/rate.dat";
 		file.open(filename.c_str());
        if (!file) {std::cerr << "Error! Unable to write " << filename << "\n"; exit(1);}
        file << _rate << "\n";
        file.close();
    }
	
	void ForwardFlux::WriteFFSConfiguration(Snapshot* snapshot, FFSConfigID& ffsconfig, bool wassuccess)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();
        //unsigned natoms = snapshot->GetNumAtoms();
		//const auto& dumpfilename = snapshot->GetSnapshotID();

        // Write the dump file out
		std::ofstream file;
        std::string filename;
        if (wassuccess)
          filename = _output_directory + "/l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + ".dat";
        else
          filename = _output_directory + "/fail-l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + ".dat";

 		file.open(filename.c_str());
        if(!file)
		{
			std::cerr << "Error! Unable to write " << filename << "\n"; 
			MPI_Abort(world_, EXIT_FAILURE);
		}

        //first line gives ID of where it came from
        file << ffsconfig.lprev << " " << ffsconfig.nprev << " " << ffsconfig.aprev << "\n";

        // Then write positions and velocities
 		for(size_t i = 0; i< atomID.size(); i++)
 		//for(size_t i = 0; i< natoms; i++)
 		{
 			file<<atomID[i]<<" ";
 			file<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" ";
 			file<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
		}
    }

    void ForwardFlux::OpenTrajectoryFile(std::ofstream& file)
	{
        FFSConfigID ffsconfig = myFFSConfigID;
        std::string filename = _output_directory + "/traj-l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + "-a" + std::to_string(ffsconfig.a) + ".xyz";
 		file.open(filename.c_str());
        if (!file)
		{
			std::cerr << "Error! Unable to write " << filename << "\n"; 
			MPI_Abort(world_, EXIT_FAILURE);
		}
    }

    void ForwardFlux::AppendTrajectoryFile(Snapshot* snapshot, std::ofstream& file)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();

        //first line gives number of atoms
        file << atomID.size() << "\n\n";

        // Then write positions and velocities
 		for(size_t i = 0; i< atomID.size(); i++)
 		{
 			file<<atomID[i]<<" ";
 			file<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" ";
 			file<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
		}
    }

	void ForwardFlux::ReadFFSConfiguration(Snapshot* snapshot, FFSConfigID& ffsconfig, bool wassuccess)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		//auto& ID = snapshot->GetSnapshotID();

        std::ifstream file; 
        std::string filename;
        if (wassuccess)
          filename =  _output_directory + "/l"+ std::to_string(ffsconfig.l) +"-n"+ std::to_string(ffsconfig.n) + ".dat";
        else
          filename =  _output_directory + "/fail-l"+ std::to_string(ffsconfig.l) +"-n"+ std::to_string(ffsconfig.n) + ".dat";

        std::string line;
        file.open(filename);
        if(!file) 
		{
			std::cerr << "Error! Unable to read " << filename << "\n"; 
			MPI_Abort(world_, EXIT_FAILURE);
		}

        unsigned int line_count = 0; 
        while(!std::getline(file,line).eof())
		{
            int atomindex = -1;
            
            //parse line into tokens
			std::string buf; // Have a buffer string
			std::stringstream ss(line); // Insert the string into a stream
			std::vector<std::string> tokens; // Create vector to hold our words
			while (ss >> buf)
			    tokens.push_back(buf);

            // first line contains the previous ffsconfig information
            if ((line_count == 0) && (tokens.size() == 3))
			{
                ffsconfig.lprev = std::stoi(tokens[0]);
                ffsconfig.nprev = std::stoi(tokens[1]);
                ffsconfig.aprev = std::stoi(tokens[2]);
            }
            // all other lines contain config information
            else if ((line_count != 0) && (tokens.size() == 7))
			{
                //FIXME: try using snapshot->GetLocalIndex()

                //copied from Ben's previous implementation
                for(size_t i=0; i < atomID.size(); i++)
                {
                    if(atomID[i] == std::stoi(tokens[0]))
                        atomindex = i;
                }

                if(atomindex < 0)
                {
                    std::cout<<"error, could not locate atomID "<<tokens[0]<<" from dumpfile"<<std::endl;
					MPI_Abort(world_, EXIT_FAILURE);
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
            // else throw error
            else
			{
				std::cout<<"ERROR: incorrect line format in "<< filename <<" on line" << line_count << ":\n";
				std::cout<<line<<std::endl;
				MPI_Abort(world_, EXIT_FAILURE);
            }
            line_count++;
        }
        file.close();
	}

    void ForwardFlux::PrintQueue()
	{
        for (unsigned int i =0 ;i < FFSConfigIDQueue.size(); i++)
		{
            std::cout << i <<" "
                      <<FFSConfigIDQueue[i].l <<" "
                      <<FFSConfigIDQueue[i].n <<" "
                      <<FFSConfigIDQueue[i].a <<" "
                      <<FFSConfigIDQueue[i].lprev <<" "
                      <<FFSConfigIDQueue[i].nprev <<" "
                      <<FFSConfigIDQueue[i].aprev <<"\n";
        }
    }

    void ForwardFlux::PopQueueMPI(Snapshot* snapshot, const CVList& cvs, unsigned int shouldpop_local)
	{
        std::vector<unsigned int> shouldpop (world_.size(),0);
        MPI_Allgather(&shouldpop_local,1,MPI_UNSIGNED,shouldpop.data(),1,MPI_UNSIGNED,world_);

        // I don't pass the queue information between procs, but I do syncronize 'shouldpop'.
        // As a result, all proc should have the same queue throughout the simulation.
		for (int i=0;i<world_.size();i++)
		{
        	if (shouldpop[i] == true)
		  	{ 
				if (i == world_.rank())
				{ //if rank matches read and pop
					if (!FFSConfigIDQueue.empty())
					{ //if queue has tasks
						myFFSConfigID = FFSConfigIDQueue.front();
						ReadFFSConfiguration (snapshot, myFFSConfigID,true);
						
						//open new trajectory file, write first frame
						if (_saveTrajectories)
						{ 
							OpenTrajectoryFile(_trajectory_file);
							AppendTrajectoryFile(snapshot,_trajectory_file);
						}

						//Trigger a rebuild of the CVs since we reset the positions
						cvs[0]->Evaluate(*snapshot);
						_cvvalue = cvs[0]->GetValue();
						_cvvalue_previous = _cvvalue;

						_pop_tried_but_empty_queue = false;
						FFSConfigIDQueue.pop_front();
					}
					else
					{ //queue is empty, need to wait for new tasks to come in
						_pop_tried_but_empty_queue = true;
					}
				}
				else
				{ //else if rank doesnt match, just pop 
					if (!FFSConfigIDQueue.empty())
						FFSConfigIDQueue.pop_front();
				}
			}
        }
        // ==============================
    }

    void ForwardFlux::ComputeCommittorProbability(Snapshot *snapshot)
	{
        // two dim vectors, first dim is lambda, second is N
        // this counts the number of atempts from this config eventually reached A (for nA) or B (for nB)
        // used to compute _pB via:  _pB = nB / (nA + nB)
        std::vector<std::vector<unsigned int>> nA;
        std::vector<std::vector<unsigned int>> nB;
        nA.resize(_ninterfaces);
        nB.resize(_ninterfaces);
        for (unsigned int i=0; i<_ninterfaces;i++)
		{
            nA[i].resize(_N[i],0);
            nB[i].resize(_N[i],0);
        }


        //Snapshot* snapshot_empty;
        FFSConfigID ffsconfig;
        
        //populate nB
        unsigned int nsuccess = _N[_ninterfaces-1]; //note -1
        for(unsigned int i = 0; i < nsuccess ; i++)
		{
        	ffsconfig.l = _ninterfaces-1;
			ffsconfig.n = i;
			ffsconfig.a = 0;
			bool flag = true;
			//recursively trace successful path and update all the configs it came from
			while (flag)
			{
				if(ffsconfig.l == 0) flag = false;

				//this is just to populate ffsconfig.{lprev,nprev,aprev}
				ReadFFSConfiguration(snapshot,ffsconfig,true);

				//update nB 
				nB[ffsconfig.l][ffsconfig.n]++;
				ffsconfig.l = ffsconfig.lprev;
				ffsconfig.n = ffsconfig.nprev;
				ffsconfig.a = ffsconfig.aprev;
			}
        }

        //now populate nA (very similar to nB)
        unsigned int nfail = _nfailure_total; 
        for(unsigned int i = 0; i < nfail ; i++)
		{
        	ffsconfig.l = 0; //only fail at lambda0 in absence of pruning
			ffsconfig.n = i;
			ffsconfig.a = 0;
			bool flag = true;
			//this is just to populate ffsconfig.{lprev,nprev,aprev}
			ReadFFSConfiguration(snapshot,ffsconfig,false);

			//recursively trace successful path and update all the configs it came from
			while (flag)
			{
				if ((ffsconfig.l == 0) && (ffsconfig.lprev==0))
					flag = false;

				//update nA 
				nA[ffsconfig.l][ffsconfig.n]++;

				ReadFFSConfiguration(snapshot,ffsconfig,true);
				ffsconfig.l = ffsconfig.lprev;
				ffsconfig.n = ffsconfig.nprev;
				ffsconfig.a = ffsconfig.aprev;
			}
        }

        //Compute _pB
        _pB.resize(_ninterfaces);
        unsigned int Nmax = 0;
        for(unsigned int i = 0; i<_ninterfaces ; i++)
		{
        	_pB[i].resize(_N[i],0);
          	for(unsigned int j = 0; j<_N[i] ; j++)
			{
            	int denom = nA[i][j] + nB[i][j];
            	if(denom != 0)
              		_pB[i][j] = (double)nB[i][j] / denom;
          	}
          	if (_N[i] > Nmax)
			  Nmax = _N[i];
        }

        //print pB
        //rows are lambda, cols are n. Thats why its looks complicated
        std::ofstream file;
        std::string filename;
        filename = _output_directory +"/commitor_probabilities.dat";
        file.open(filename.c_str());
        if (!file)
		{
			std::cerr << "Error! Unable to write " << filename << "\n"; 
			MPI_Abort(world_, EXIT_FAILURE);
		}

        for(unsigned int i = 0; i<Nmax ; i++)
		{
        	for(unsigned int j = 0; j<_ninterfaces ; j++)
			{
				if (i < _N[j])
					file << _pB[j][i] << " ";
				else
					file << "none ";
          	}
          	file << "\n";
		}
    }
    
    void ForwardFlux::ReconstructTrajectories(Snapshot *snapshot)
	{
        //this only reconstructs successful trajectories 
        //its pretty straightforward to reconstruct failed trajectories, but I dont do it here
        int nsuccess = _S[_ninterfaces-2]; //note -2, no attempts from _ninterfaces-1
        
        //Snapshot* snapshot_empty;
        FFSConfigID ffsconfig;

        for(int i = 0; i < nsuccess ; i++)
		{
			std::deque<FFSConfigID> path;

			ffsconfig.l = _ninterfaces-1;
			ffsconfig.n = i;
			ffsconfig.a = 0;
			bool flag = true;

			while (flag)
			{
				if (ffsconfig.l == 0){ flag = false;}

				//this is just to populate ffsconfig.{lprev,nprev,aprev}
				//ReadFFSConfiguration(snapshot_empty,ffsconfig);
				ReadFFSConfiguration(snapshot,ffsconfig,true);
				
				//path.emplace_front(ffsconfig); //not sure if new constructor will work
				path.emplace_front(ffsconfig.l,ffsconfig.n, ffsconfig.a, ffsconfig.lprev, ffsconfig.nprev, ffsconfig.aprev); 

				ffsconfig.l = ffsconfig.lprev;
				ffsconfig.n = ffsconfig.nprev;
				ffsconfig.a = ffsconfig.aprev;
			}
			//the last element is the last lambda, which doesn't have a trajectory, so delete it
			path.pop_back();

			//now path should contain all of the FFSConfigID's from B back to A
			//reverse pop it and splice all traj- files into a new traj-full- file

			//output file
			std::ofstream ofile;
			std::string ofilename;
			ofilename = _output_directory +"/traj-full-" + std::to_string(i) + ".xyz";
			ofile.open(ofilename.c_str());
			if(!ofile)
			{
				std::cerr << "Error! Unable to write " << ofilename << "\n"; 
				MPI_Abort(world_, EXIT_FAILURE);
			}

			while(!path.empty())
			{
				ffsconfig = path.front();
				path.pop_front();

				//input file
				std::ifstream ifile;
				std::string ifilename;

				ifilename = _output_directory + "/traj-l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + "-a" + std::to_string(ffsconfig.a) + ".xyz";

				ifile.open(ifilename.c_str());
				if(!ifile)
				{
					std::cerr << "Error! Unable to read " << ifilename << "\n";
					MPI_Abort(world_, EXIT_FAILURE);
				}
				//write entire ifile to ofile
				std::string line;
				while(!std::getline(ifile,line).eof())
					ofile << line << "\n";
			}
			ofile.close();
		}
	}	

	ForwardFlux* ForwardFlux::Build(const Json::Value& json, 
	                                const MPI_Comm& world,
	                                const MPI_Comm& comm,
	                                const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		CharReaderBuilder rbuilder;
		CharReader* reader = rbuilder.newCharReader();

		reader->parse(JsonSchema::ForwardFluxMethod.c_str(),
		              JsonSchema::ForwardFluxMethod.c_str() + JsonSchema::ForwardFluxMethod.size(),
		              &schema, nullptr);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());
		
		double ninterfaces = json.get("nInterfaces", 2).asDouble();
		std::vector<double> interfaces;
		for(auto& s : json["interfaces"])
			interfaces.push_back(s.asDouble());

		std::vector<unsigned int> M;
		for(auto& s : json["trials"])
			M.push_back(s.asInt());
		
		if ((ninterfaces != interfaces.size()) || (ninterfaces != M.size()))
			throw BuildException({"The size of \"interfaces\" and \"trials\" must be equal to \"nInterfaces\". See documentation for more information"});

		auto N0Target = json.get("N0Target", 1).asInt();
		auto NLastSuccessful = json.get("NLastSuccessful", 1).asInt();
		auto initialFluxFlag = json.get("computeInitialFlux", true).asBool();
		auto saveTrajectories = json.get("saveTrajectories", true).asBool();
		auto currentInterface = json.get("currentInterface", 0).asInt();
		auto freq = json.get("frequency", 1).asInt();
		auto flavor = json.get("flavor", "none").asString();
		auto output_directory = json.get("outputDirectoryName", "FFSoutput").asString();

		// TODO: Implement multiple processes per walker
		if(mxx::comm(comm).size() > 1)
		{
			throw BuildException({"Forward Flux currently only works with 1 process per walker."});
		}

		if(flavor == "DirectForwardFlux")
		{
			return new DirectForwardFlux(world, comm, ninterfaces, interfaces, N0Target, NLastSuccessful, M, initialFluxFlag, saveTrajectories, currentInterface, output_directory, freq);            	
		}
		else 
		{
			throw BuildException({"Unknown flavor of forward flux. The options are \"DirectForwardFlux\""});
		}
	}
}
