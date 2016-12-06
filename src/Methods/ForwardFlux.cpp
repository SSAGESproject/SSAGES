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
#include <queue>

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
        if (_iteration % 1 != 0) return;

        // if _computefluxA0
        if (_initialFluxFlag){
            ComputeInitialFlux(snapshot,cvs); 
            InitializeQueue(snapshot,cvs);
        }
        // Else check the FFS interfaces
        else{
            CheckForInterfaceCrossings(snapshot,cvs);
          
        }
        // Other modes?

    }

	void ForwardFlux::PostSimulation(Snapshot* snapshot, const CVList& cvs)
	{
		
	}
    
    int ForwardFlux::HasCrossedInterface(double current, double prev, unsigned int i){
        double interface_location = _interfaces[i];
        if ((prev < interface_location) && (current >= interface_location))
            return 1;
        else if ((prev >= interface_location) && (current < interface_location))
            return -1;
        else
            return 0;
    }

    bool ForwardFlux::HasReturnedToA(double current){
        double interface_location = _interfaces[0];
        if (current < interface_location) return true;
        else return false;
    }

	void ForwardFlux::ComputeInitialFlux(Snapshot* snapshot, const CVList& cvs){
        //check if I've crossed the first interface (lambda 0)

        //need to sync variables between processors
    }

    void ForwardFlux::CheckForInterfaceCrossings(Snapshot* snapshot, const CVList& cvs){

        //QUESTION: Whats the difference between _world and _comm?
        //For now I'll use _world for everything. But if each driver uses multiple procs then I suspect that this will be wrong.

        _cvvalue = cvs[0]->GetValue();
        //check if I've crossed the next interface
        bool hasreturned = HasReturnedToA(_cvvalue);
        int hascrossed = HasCrossedInterface(_cvvalue, _cvvalue_previous, _current_interface + 1);
        bool fail_local=false,success_local=false;
        bool *successes = new bool (_world.size());
        bool *failures = new bool (_world.size());
        //bool vectors in MPI were strange
        //std::vector<bool> successes;
        //std::vector<bool> failures;
        //successes.resize(_world.size());
        //failures.resize(_world.size());
        
        if (!_succeeded_but_empty_queue && !_failed_but_empty_queue){ 
            // make sure this isnt a zombie trajectory that previously failed or succeeded and is just waiting for the queue to get more jobs
            if (hasreturned){ 
              fail_local=true;
            }
            else if (hascrossed == 1){
              success_local=true;
            }
            else if (hascrossed == -1){
              //this should never happen if the interfaces are non-intersecting, it would be wise to throw an error here though
            }
            else{
              //not sure if anything should needs to be done here        
            }
        }

        //for each traj that crossed to lambda+1 need to write it to disk (FFSConfigurationFile)
        //MPIAllgather success_local into successes
        MPI_Allgather(&success_local,1,MPI::BOOL,successes,1,MPI::BOOL,_world);
        MPI_Allgather(&fail_local,1,MPI::BOOL,failures,1,MPI::BOOL,_world);
       
        int success_count = 0, fail_count = 0;
        // I dont pass the queue information between procs but I do syncronize 'successes' and 'failures'
        //   as a reuslt all proc should have the same queue throughout the simulation
        for (int i=0;i<_world.size();i++){
          if (successes[i] == true){ 
            if (i == _world.rank()){
              // write config to lambda+1
              int l,n,a,lprev,nprev,aprev;
              lprev = myFFSConfigID.l;
              nprev = myFFSConfigID.n;
              aprev = myFFSConfigID.a;
              //update ffsconfigid's l,n,a
              l = lprev + 1;
              n = _S[_current_interface] + 1 + success_count;
              a = 0; //in DFFS, everyone gets one attempt (unless randomly you choose the same config to shoot from twice...I should look into whether this is allowed). At the very least however, a=0 to start with the possibility that it will be >0 if same config is chosen twice.
              
              FFSConfigID newid = FFSConfigID(l,n,a,lprev,nprev,aprev);
              WriteFFSConfiguration(snapshot,newid);
            }
            success_count++;
          }
          if (failures[i] == true){ 
            fail_count++;
          }
         
        }
        
        //update the number of successes and attempts, same for all proc since Allgathered 'successes' and 'failures'
        _S[_current_interface] += success_count;
        _attempts[_current_interface] += success_count + fail_count;
        // ^ I dont like storing attempts this way (as is its only when they finish). _attempts should also include jobs in progress (i.e. jobs currently running). THINK ABOUT THIS!.
        
        // Check if this interface is finished, if so add new tasks to queue
        if (_S[_current_interface] == _M[_current_interface]) {

          //for DFFS
          unsigned int npicks = _M[_current_interface+1];
          std::vector<unsigned int> picks;
          picks.resize(npicks);

          if (_world.rank() == 0){
            std::uniform_int_distribution<int> distribution(0,_S[_current_interface]-1);
            for (int i=0; i < npicks ; i++){
               picks[i] = distribution(_generator);
            }
          }
          MPI_Bcast(picks.data(),npicks,MPI::UNSIGNED,0,_world);


          //each proc adds to the queue
          for (int i=0; i < npicks ; i++){
            int l,n,a,lprev,nprev,aprev;
            lprev = myFFSConfigID.l;
            nprev = myFFSConfigID.n;
            aprev = myFFSConfigID.a;
            //update ffsconfigid's l,n,a
            l = lprev + 1;
            n = i;
            a = 0; //dont correct for choosing the same config twice (should I?)
            FFSConfigIDQueue.emplace(l,n,a,lprev,nprev,aprev);
          }
        }

        // if succeeded or failed (or zombie job), get a new config from the queue...but need to be careful that no two procs get the same config

        // Need to account for zombie jobs that are waiting for a new config
        bool shouldpop_local = false;        
        bool *shouldpop = new bool(_world.size());
        //std::vector<bool> shouldpop;
        //shouldpop.resize(_world.size());
        if (success_local || fail_local || _succeeded_but_empty_queue || _failed_but_empty_queue){
          shouldpop_local = true;
        }
        MPI_Allgather(&shouldpop_local,1,MPI::BOOL,shouldpop,1,MPI::BOOL,_world);

        // I dont pass the queue information between procs but I do syncronize 'shouldpop'
        //   as a reuslt all proc should have the same queue throughout the simulation
        for (int i=0;i<_world.size();i++){
          if (shouldpop[i] == true){ 
            if (i == _world.rank()){ //if rank matches read and pop
              if (!FFSConfigIDQueue.empty()){ //if queue has tasks
                   myFFSConfigID = FFSConfigIDQueue.front();
                   ReadFFSConfiguration (snapshot, myFFSConfigID);
                   _succeeded_but_empty_queue = false;
                   _failed_but_empty_queue = false;
                  FFSConfigIDQueue.pop();
              }
              else{ //queue is empty, need to wait for new tasks to come in
                  if (successes[i] == true)     _succeeded_but_empty_queue = true;
                  else if (failures[i] == true) _failed_but_empty_queue = true;
              }
            }
            else{ //else if rank doesnt match, just pop 
              if (!FFSConfigIDQueue.empty()){
                  FFSConfigIDQueue.pop();
              }
            }
          }
        }

        
        //Anything else to update across mpi?


        _cvvalue_previous = _cvvalue;

        //delete[] successes;
        //delete[] failures,shouldpop;
           
    }

	void ForwardFlux::ComputeTransitionProbabilities(){

    }

	
	void ForwardFlux::WriteFFSConfiguration(Snapshot* snapshot, FFSConfigID ffsconfig)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();
		const auto& dumpfilename = snapshot->GetSnapshotID();

    }


	void ForwardFlux::ReadFFSConfiguration(Snapshot* snapshot, FFSConfigID ffsconfig)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		auto& ID = snapshot->GetSnapshotID();

	}

    void ForwardFlux::InitializeQueue(Snapshot* snapshot, const CVList &cvs){

        unsigned int npicks = _M[0];
        std::vector<unsigned int> picks;
        picks.resize(npicks);

        if (_world.rank() == 0){
          std::uniform_int_distribution<int> distribution(0,_N0-1);
          for (int i=0; i < npicks ; i++){
             picks[i] = distribution(_generator);
          }
        }
        MPI_Bcast(picks.data(),npicks,MPI::UNSIGNED,0,_world);


        //each proc adds to the queue
        for (int i=0; i < npicks ; i++){
          int l,n,a,lprev,nprev,aprev;
          FFSConfigID *myconfig = &Lambda0ConfigLibrary[picks[i]];
          lprev = myconfig->l;
          nprev = myconfig->n;
          aprev = myconfig->a;
          //update ffsconfigid's l,n,a
          // current = previous, thats how you know you're lambda0
          l = lprev;
          n = nprev;
          a = aprev; 
          FFSConfigIDQueue.emplace(l,n,a,lprev,nprev,aprev);
        }         

    }
    

	
}

