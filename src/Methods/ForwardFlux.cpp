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
        if (_initialFluxFlag){
            ComputeInitialFlux(snapshot,cvs); 
        }
        // Else normal forward flux
        else{
            CheckForInterfaceCrossings(snapshot,cvs);
          
        }
        // Other modes?

    }

	void ForwardFlux::PostSimulation(Snapshot*, const CVList&)
	{
		
	}
    
    int HasCrossedInterface(double current, double prev, unsigned int i){
        double interface_location = _interfaces[i]
        if ((prev < interface_location) && (current >= interface_location))
            return 1;
        else if ((prev >= interface_location) && (current < interface_location))
            return -1;
        else
            return 0;
    }

    bool ForwardFlux::HasReturnedToA(double current){
        double interface_location = _interfaces[0]
        if (current < interface_location) return true;
        else return false;
    }

	void ForwardFlux::ComputeInitialFlux(){
        //check if I've crossed the first interface (lambda 0)

        //need to sync variables between processors
    }

    void ForwardFlux::CheckForInterfaceCrossings(Snapshot* snapshot, CVList& cvs){

        _cvvalue = cvs[0].GetValue()
        //check if I've crossed the next interface
        bool hasreturend = HasReturnedToA(_cvvalue);
        int hascrossed = HasCrossedInterface(_cvvalue, _cv_value_previous, _currentinterface + 1);
        int nfail_local=0,nsuccess_local=0;
        //std::vector<int> nfail,nsuccess;
        //nfail.resize(SIZEWORLD)
        //nsuccess.resize(SIZEWORLD)
        
        bool shouldIpop_local = false;
        std:vector<bool> shouldIpop;
        shouldIpop.resize(SIZEWORLD);
        std::vector<int> poporder; //order of processors to pop off queue
        int myplaceinline; //where current processor is in poporder

        if (hasreturned){
          shouldIpop = true;
          fail_local=1;
        }
        else if (hascrossed == 1){
          shouldIpop = true;
          success_local=1;
        }
        else if (hascrossed == -1){
          //this should never happen if the interfaces are non-intersecting
        }
        else{
        }

        //for each traj that crossed to lambda+1 need to write it to disk (FFSConfigurationFile)
        std::vector<bool> successes;
        MPIGather success_local into successes

        for (i=0;i<SIZEWORLD;i++){
          if ((shouldIpop[i] == true){ 
            if (i == MYRANK)){
              // write config to lambda+1
              int l,n,a,l_prev,n_prev,a_prev;
              l = myFFSConfigID.getl()
              n = myFFSConfigID.getn()
              a = myFFSConfigID.geta()
              FFSConfigID newid = FFSConfigID(l,n,a,
              UPDATE FFS Config ID to lambda+1
              WriteFFSConfiguration(snapshot,FFSConfigID);
            }
          }
        }
        
        MPI_Reduce _S[_current_interface]
        //replace if statement with a generalized function like "takeAction" that
        //  - given the number of successes see if more simulations need to be spawned
        //or perhaps post_integration will just have to be different for each FFS flavor
        if (_S[_current_interface] == _M[_current_interface]) {
          AddNewIDsToQueue();
          //FFSConfigIDQueue.emplace(NEW AND CORRECT ID INFO);
        }

        // if returned or crossed, get a new config from the queue. but need to be careful that no two procs get the same config
        MPI_Gather into shouldIpop
         
        for (i=0;i<SIZEWORLD;i++){
          if ((shouldIpop[i] == true){ 
              if (i == MYRANK)){
                  myFFSConfigID = FFSConfigIDQueue.front()
              }
              FFSConfigIDQueue.pop()
          }
        }

        MPI_Gather into nsuccess
        MPI_Gather into nfail

        //MPI_Reduce _S

          //update snapshot with new traj
          myFFSConfigID = FFSConfigIDQueue.pop();
        
        //if so update some relevant quantities and mpi them across procs

        //assign myFFSConfigID = NULL

        // reassign configs that reached interface a new FFSConfigID

        _cvvalue_previous = _cvvalue;
           
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
    void ForwardFlux::AddNewIDsToQueue(){
        //for DFFS
        if (rank == 0){
          roll M[nextinterface] random numbers (these are the simulations at lambda_i+1 to pick)
        }
        MPI_Broadcast these random numbers

        each proc adds to the queue


        //for CBGFFS
        similar to DFFS?

        //for BGFFS
        will need to have broadcast ID information of the successful config that is spawning new configs



    }

	
}

