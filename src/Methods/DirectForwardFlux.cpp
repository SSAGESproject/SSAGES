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
#include "DirectForwardFlux.h"
#include "FileContents.h"
#include "CVs/CVManager.h"
#include <iostream>
#include <random>
#include <queue>

namespace SSAGES
{
	void DirectForwardFlux::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		//check if we want to check FFS interfaces this timestep
		//for now, do it every time step
		if (iteration_ % 1 != 0) return;
		
		auto cvs = cvmanager.GetCVs(cvmask_);
		// check the structure at the beginning of the simulation
		if (iteration_ == 0)
			CheckInitialStructure(cvs);

		// if _computefluxA0
		if (_initialFluxFlag)
		{
			ComputeInitialFlux(snapshot,cvs); 
			if (!_initialFluxFlag)
			{ //only enter here once
				InitializeQueue(snapshot,cvs);
			  	PrintQueue();
			}
		}
		else if (initializeQueueFlag)
		{
			InitializeQueue(snapshot,cvs);
			PrintQueue();
		}
		// Else check the FFS interfaces
		else
		{
			CheckForInterfaceCrossings(snapshot, cvmanager); 
			//FluxBruteForce(snapshot,cvs);
		}
		// Other modes?
	}

	void DirectForwardFlux::CheckForInterfaceCrossings(Snapshot* snapshot, const CVManager& cvmanager)
	{
		//This is the main FFS method. The magic happens here!
		auto cvs = cvmanager.GetCVs(cvmask_);
		
		//QUESTION: Whats the difference between world_ and _comm?
		//For now I'll use world_ for everything. But if each driver uses multiple procs then I suspect that this will be wrong.
		_cvvalue = cvs[0]->GetValue();
		
		//// for debugging
		//if ((myFFSConfigID.l == 2) && (myFFSConfigID.n == 20) && (myFFSConfigID.a ==1)){
		//    std::cout << "stop here\n";
		//}
		//if (myFFSConfigID.l != _current_interface){
		//    std::cout << "stop here\n";
		//}

		//check if I've crossed the next interface
		bool hasreturned = HasReturnedToA(_cvvalue);
		int hascrossed = HasCrossedInterface(_cvvalue, _cvvalue_previous, _current_interface + 1);
		unsigned int fail_local=false,success_local=false;
		//std::vector<bool> in MPI were strange, ended up using std::vector<unsigned int>
		std::vector<unsigned int> successes (world_.size());
		std::vector<unsigned int> failures (world_.size());

		if (!_pop_tried_but_empty_queue)
		{ 
			// make sure this isnt a zombie trajectory that previously failed or succeeded and is just waiting for the queue to get more jobs
			if (hasreturned)
				fail_local=true;
			else if (hascrossed == 1)
				success_local=true;
			else if (hascrossed == -1)
			{
			  //this should never happen if the interfaces are non-intersecting, it would be wise to throw an error here though
			  std::cerr << "Trajectory l"<<myFFSConfigID.l<<"-n"<<myFFSConfigID.n<<"-a"<<myFFSConfigID.a<<" crossed backwards! This shouldnt happen!" << std::endl;
			  //world_.abort(EXIT_FAILURE);
			}
			else
			{
			  //not sure if anything should needs to be done here        
			}
		}

		//for each traj that crossed to lambda+1 need to write it to disk (FFSConfigurationFile)
		//MPIAllgather success_local into successes
		MPI_Allgather(&success_local,1,MPI_UNSIGNED,successes.data(),1,MPI_UNSIGNED,world_);
		MPI_Allgather(&fail_local,1,MPI_UNSIGNED,failures.data(),1,MPI_UNSIGNED,world_);
	   
		int success_count = 0, fail_count = 0;
		// I dont pass the queue information between procs but I do syncronize 'successes' and 'failures'
		//   as a reuslt all proc should have the same queue throughout the simulation
		for (int i=0;i<world_.size();i++)
		{
			int l,n,a,lprev,nprev,aprev;
			// write config to lambda+1
			lprev = myFFSConfigID.l;
			nprev = myFFSConfigID.n;
			aprev = myFFSConfigID.a;

			if (successes[i] == true)
			{ 
				if (i == world_.rank())
				{
					//update ffsconfigid's l,n,a
					l = _current_interface + 1;
					n = _S[_current_interface] + success_count;
					a = 0;
					FFSConfigID newid = FFSConfigID(l,n,a,lprev,nprev,aprev);
					WriteFFSConfiguration(snapshot,newid,1);
				}
				success_count++;
			}
			if (failures[i] == true)
			{ 
				if (i == world_.rank())
				{
					//update ffsconfigid's l,n,a
					l = 0; //only fail at lambda 0, 
					n = _nfailure_total + fail_count;
					a = 0;
					FFSConfigID newid = FFSConfigID(l,n,a,lprev,nprev,aprev);
					WriteFFSConfiguration(snapshot,newid,0);
				}
				fail_count++;
			}
		}

		//update trajectories
		if (_saveTrajectories)
		{
		  	if (!_pop_tried_but_empty_queue) //dont update trajectory if zombie job
				AppendTrajectoryFile(snapshot,_trajectory_file);
		  	if (success_local || fail_local) //if finished then close it
				_trajectory_file.close();
		}
		
		//update the number of successes and attempts, same for all proc since Allgathered 'successes' and 'failures'
		_S[_current_interface] += success_count;
		_A[_current_interface] += success_count + fail_count;
		// ^ I dont like storing attempts this way (as is its only when they finish). _A should also include jobs in progress (i.e. jobs currently running). THINK ABOUT THIS!.
		_nfailure_total += fail_count;

		//------------------------------
		//print info
		if ((success_local) || (fail_local))
		{
			std::cout << "Iteration: "<< iteration_ << ", proc " << world_.rank() << "\n";
		  	if (success_local)
			{
				std::cout << "Successful attempt from interface " << _current_interface 
						<< " l"<<myFFSConfigID.l<<"-n"<<myFFSConfigID.n<<"-a"<<myFFSConfigID.a
						<< " (cvvalue_previous: " << _cvvalue_previous << " cvvalue " << _cvvalue << " interface["<<_current_interface+1<<"] = " << _interfaces[_current_interface+1] << "\n";
		  	}
		  	if (fail_local)
			{
				std::cout << "Failed attempt from interface " << _current_interface 
						<< " l"<<myFFSConfigID.l<<"-n"<<myFFSConfigID.n<<"-a"<<myFFSConfigID.a
						<< " (cvvalue_previous: " << _cvvalue_previous << " cvvalue " << _cvvalue << " interface[0] = "<< _interfaces[0] << "\n"
						<< "nfailuretotal: " << _nfailure_total << "\n";
		  	}
			
			std::cout << "A: ";
			for (auto a : _A)
				std::cout << a << " "; 
			std::cout << "\n";

			std::cout << "S: ";
			for (auto s : _S)
				std::cout << s << " ";
			std::cout << "\n";
			
			std::cout << "M: ";
			for (auto m : _M)
				std::cout << m << " ";
			std::cout << "\n";
		}
		//------------------------------
	   
		//create new funciton here? (call it SetupNextInterface()
		// Check if this interface is finished, if so add new tasks to queue, and increment _current_interface
		if (_A[_current_interface] >= _M[_current_interface])
		{
			//set N and P
			_N[_current_interface+1] = _S[_current_interface];
			_P[_current_interface] = (double) _S[_current_interface] / _A[_current_interface];

			if (_current_interface+2 < _ninterfaces)
			{
				_current_interface += 1;
				//_N[_current_interface] = _S[_current_interface-1];

				if (_N[_current_interface] == 0)
				{
					std::cerr << "Error! No successes from interface " << _current_interface-1 << " to " << _current_interface <<"! Try turning up M["<<_current_interface-1<<"] or spacing interfaces closer together.\n";
					MPI_Abort(world_, EXIT_FAILURE);
				}

				//for DFFS
				unsigned int npicks = _M[_current_interface];
				std::vector<unsigned int> picks;
				picks.resize(npicks);
				if (IsMasterRank(world_))
				{
					std::uniform_int_distribution<int> distribution(0,_N[_current_interface]-1);
					for (unsigned int i=0; i < npicks ; i++)
						picks[i] = distribution(_generator);
				}
				MPI_Bcast(picks.data(),npicks,MPI_UNSIGNED,0,world_);

				//each proc adds to the queue
				//set correct attempt index if a given ID is picked twice
				std::vector<unsigned int> attempt_count;
				attempt_count.resize(_N[_current_interface],0);

				for (unsigned int i=0; i < npicks ; i++)
				{
					unsigned int mypick = picks[i];
					int l,n,a,lprev,nprev,aprev;
					lprev = myFFSConfigID.l;
					nprev = myFFSConfigID.n;
					aprev = myFFSConfigID.a;
					//update ffsconfigid's l,n,a
					l = _current_interface;
					n = mypick;
					a = attempt_count[mypick]; 
					attempt_count[mypick]++; //this updates attempt number if same config is picked twice

					FFSConfigIDQueue.emplace_back(l,n,a,lprev,nprev,aprev);
				}
			}
			else
			{
				std::cout << "DFFS should be finished here, do something special? like exit?\n";
				//Hythem said this is "acceptable" until the code is changed
				PostSimulation(snapshot, cvmanager);
			}
		}

		// if succeeded or failed (or zombie job), get a new config from the queue...but need to be careful that no two procs get the same config

		// Need to account for zombie jobs that are waiting for a new config
		unsigned int shouldpop_local = false;        
		std::vector<unsigned int> shouldpop(world_.size(),0);

		if (success_local || fail_local || _pop_tried_but_empty_queue)
		  shouldpop_local = true;
		
		//Pop the queue
		// Need to perform mpi call so that all proc pop the queue in the same way
		PopQueueMPI(snapshot,cvs, shouldpop_local);
				
		//Anything else to update across mpi?
		//clean up
		_cvvalue_previous = _cvvalue;
		iteration_++;
	}

	void DirectForwardFlux::InitializeQueue(Snapshot* snapshot, const CVList &cvs)
	{
		unsigned int npicks = _M[0];
		std::vector<unsigned int> picks;
		picks.resize(npicks);
		
		// if the initial flux has already been realized or
		// we've got more than _N0Target configurations for zero surface 		
		if (_current_interface == 0)
			_N[_current_interface]=_N0Target;

		//if the calculation starts from a non-zero surface (for restart)
		if (_N[_current_interface] == 0)
			_N[_current_interface]=_NLastSuccessful;

		if(IsMasterRank(world_))
		{
			std::uniform_int_distribution<int> distribution(0,_N[_current_interface]-1);
		  	for (unsigned int i=0; i < npicks ; i++)
				picks[i] = distribution(_generator);
		}
		MPI_Bcast(picks.data(),npicks,MPI_UNSIGNED,0,world_);

		//set correct attempt index if a given ID is picked twice
		std::vector<unsigned int> attempt_count;
		attempt_count.resize(_N[_current_interface],0);

		//each proc adds to the queue
		for (unsigned int i=0; i < npicks ; i++)
		{
			unsigned int mypick = picks[i];
			int l,n,a,lprev,nprev,aprev;
			FFSConfigID *myconfig = &Lambda0ConfigLibrary[picks[i]];
			lprev = myconfig->l;
			nprev = myconfig->n;
			aprev = myconfig->a;
			//update ffsconfigid's l,n,a
			// current = previous, thats how you know you're lambda0
			l = lprev;
			n = nprev;
			a = attempt_count[mypick]; 
			attempt_count[mypick]++; //this updates attempt number if same config is picked twice
			FFSConfigIDQueue.emplace_back(l,n,a,lprev,nprev,aprev);
		}         
		std::cout << "FFSConfigIDQueue has " << FFSConfigIDQueue.size() << " entries upon initialization\n";
		initializeQueueFlag = false;
		
		// now that queue is populated initialize tasks for all processors
		// ==============================
		unsigned int shouldpop_local = true;
		PopQueueMPI(snapshot,cvs,shouldpop_local);
	}	
}

