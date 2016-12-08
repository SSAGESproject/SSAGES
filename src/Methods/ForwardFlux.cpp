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


        std::cout << "WARNING! MAKE SURE LAMMPS GIVES A DIFFERENT RANDOM SEED TO EACH PROCESSOR, OTHERWISE EACH FFS TRAJ WILL BE IDENTICAL!\n";

        //code for setting up simple simulation and debugging
        _ninterfaces = 4;
        int i;
        _interfaces.resize(_ninterfaces);
        _interfaces[0]=-0.8;
        _interfaces[1]=-0.5;
        _interfaces[2]= 0.0;
        _interfaces[3] =0.5;

        _current_interface = 0;

        output_directory = "FFSoutput";
        _initialFluxFlag = true;
        _M.resize(_ninterfaces);
        _A.resize(_ninterfaces);
        _P.resize(_ninterfaces);
        _S.resize(_ninterfaces);
        _N.resize(_ninterfaces);
        for(i=0;i<_ninterfaces;i++) _M[i] = 5;

        _N0 = 10;
        Lambda0ConfigLibrary.resize(_N0);
          std::normal_distribution<double> distribution(0,1);
        for (i = 0; i < _N0 ; i++){
          Lambda0ConfigLibrary[i].l = 0;
          Lambda0ConfigLibrary[i].n = i;
          Lambda0ConfigLibrary[i].a = 0;
          Lambda0ConfigLibrary[i].lprev = 0;
          Lambda0ConfigLibrary[i].nprev = i;
          Lambda0ConfigLibrary[i].aprev = 0;

          FFSConfigID ffsconfig = Lambda0ConfigLibrary[i];
          // Write the dump file out
          std::ofstream file;
          std::string filename = output_directory + "/l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + "-a" + std::to_string(ffsconfig.a) + ".dat";
          file.open(filename.c_str());

          //first line gives ID of where it came from
          file << ffsconfig.lprev << " " << ffsconfig.nprev << " " << ffsconfig.aprev << "\n";
          //write position and velocity
          file << "1 -1 0 0 " << distribution(_generator) << " "  << distribution(_generator) << " 0\n";

        }
        _pop_tried_but_empty_queue = false;


    }

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
        //check if we want to check FFS interfaces this timestep
        if (_iteration % 1 != 0) return;


        // if _computefluxA0
        if (_initialFluxFlag){
            ComputeInitialFlux(snapshot,cvs); 
            if (!_initialFluxFlag){
              InitializeQueue(snapshot,cvs);
              PrintQueue();
            }
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

        //responsible for setting _initialFluxFlag = false when finished
        _initialFluxFlag = false;
    }

    void ForwardFlux::CheckForInterfaceCrossings(Snapshot* snapshot, const CVList& cvs)
    {

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
        
        if (!_pop_tried_but_empty_queue){ 
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
        _A[_current_interface] += success_count + fail_count;
        // ^ I dont like storing attempts this way (as is its only when they finish). _A should also include jobs in progress (i.e. jobs currently running). THINK ABOUT THIS!.
        
        // Check if this interface is finished, if so add new tasks to queue, and increment _current_interface
        if (_S[_current_interface] == _M[_current_interface]){
          if (_current_interface+1 != _ninterfaces){
            _current_interface += 1;
            _N[_current_interface] = _S[_current_interface-1];

            //for DFFS
            unsigned int npicks = _M[_current_interface];
            std::vector<unsigned int> picks;
            picks.resize(npicks);

            if (_world.rank() == 0){
              std::uniform_int_distribution<int> distribution(0,_N[_current_interface]-1);
              for (int i=0; i < npicks ; i++){
                 picks[i] = distribution(_generator);
              }
            }
            MPI_Bcast(picks.data(),npicks,MPI::UNSIGNED,0,_world);


            //each proc adds to the queue

            //set correct attempt index if a given ID is picked twice
            std::vector<unsigned int> attempt_count;
            attempt_count.resize(_N[_current_interface],0);

            for (int i=0; i < npicks ; i++){
              unsigned int mypick = picks[i];
              int l,n,a,lprev,nprev,aprev;
              lprev = myFFSConfigID.l;
              nprev = myFFSConfigID.n;
              aprev = myFFSConfigID.a;
              //update ffsconfigid's l,n,a
              l = lprev + 1;
              n = i;
              a = attempt_count[mypick]; 
              attempt_count[mypick]++; //this updates attempt number if same config is picked twice

              FFSConfigIDQueue.emplace_back(l,n,a,lprev,nprev,aprev);
            }
          }
          else{
            std::cout << "DFFS should be finished here, do something special? like exit?\n";
			_world.abort(EXIT_FAILURE); //more elegant solution?
          }
        }

        // if succeeded or failed (or zombie job), get a new config from the queue...but need to be careful that no two procs get the same config

        // Need to account for zombie jobs that are waiting for a new config
        bool shouldpop_local = false;        
        bool *shouldpop = new bool(_world.size());
        //std::vector<bool> shouldpop;
        //shouldpop.resize(_world.size());
        if (success_local || fail_local || _pop_tried_but_empty_queue){
          shouldpop_local = true;
        }
        
        //Pop the queue
        // Need to perform mpi call so that all proc pop the queue in the same way
        PopQueueMPI(shouldpop_local);
                
        //Anything else to update across mpi?
        

        //print info
        std::cout << "Iteration: "<< _iteration << ", proc " << _world.rank() << "\n";
        if (success_local)
          std::cout << "Successful attempt from interface " << _current_interface-1 <<" (cvvalue_previous: " << _cvvalue_previous << "cvvalue " << _cvvalue << "interface " << _interfaces[_current_interface-1] << "\n";
        if (fail_local)
          std::cout << "Failed attempt from interface " << _current_interface-1 <<" (cvvalue_previous: " << _cvvalue_previous << "cvvalue " << _cvvalue << "interface " << _interfaces[0] << "\n";
        std::cout << "A: ";
        for (auto a : _A) std::cout << a << " "; std::cout << "\n";
        std::cout << "S: ";
        for (auto s : _S) std::cout << s << " "; std::cout << "\n";
        std::cout << "M: ";
        for (auto m : _M) std::cout << m << " "; std::cout << "\n";
        
        //clean up
        _cvvalue_previous = _cvvalue;
        _iteration++;

        delete[] successes;
        delete[] failures,shouldpop;


           
    }

	void ForwardFlux::ComputeTransitionProbabilities(){

    }

	
	void ForwardFlux::WriteFFSConfiguration(Snapshot* snapshot, FFSConfigID ffsconfig)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();
		//const auto& dumpfilename = snapshot->GetSnapshotID();

        // Write the dump file out
		std::ofstream file;
        std::string filename = output_directory + "/l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + "-a" + std::to_string(ffsconfig.a) + ".dat";
 		file.open(filename.c_str());

        //first line gives ID of where it came from
        file << ffsconfig.lprev << " " << ffsconfig.nprev << " " << ffsconfig.aprev << "\n";

        // Then write positions and velocities
 		for(size_t i = 0; i< atomID.size(); i++)
 		{
 			file<<atomID[i]<<" ";
 			file<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" ";
 			file<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
		}

    }


	void ForwardFlux::ReadFFSConfiguration(Snapshot* snapshot, FFSConfigID ffsconfig)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		//auto& ID = snapshot->GetSnapshotID();

        std::ifstream file; 
        std::string filename =  output_directory + "/l"+ std::to_string(ffsconfig.l) +"-n"+ std::to_string(ffsconfig.n) +"-a"+ std::to_string(ffsconfig.a) + ".dat";
        std::string line;

        file.open(filename);
        if (!file) {std::cerr << "Error! Unable to open " << filename << "\n"; exit(1);}

        unsigned int line_count = 0; 
        while(!std::getline(file,line).eof()){

            int atomindex = -1;
            
            //parse line into tokens
			std::string buf; // Have a buffer string
			std::stringstream ss(line); // Insert the string into a stream
			std::vector<std::string> tokens; // Create vector to hold our words
			while (ss >> buf)
			    tokens.push_back(buf);
           
            // first line contains the previous ffsconfig information
            if ((line_count == 0) && (tokens.size() == 3)){
                ffsconfig.lprev = std::stoi(tokens[0]);
                ffsconfig.nprev = std::stoi(tokens[1]);
                ffsconfig.aprev = std::stoi(tokens[2]);
            }
            // all other lines contain config information
            else if ((line_count != 0) && (tokens.size() == 7)){
            
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
            // else throw error
            else{
				std::cout<<"ERROR: incorrect line format in "<< filename <<" on line" << line_count << ":\n";
				std::cout<<line<<std::endl;
				_world.abort(-1);	
            }
            line_count++;
        }
        file.close();
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
          FFSConfigIDQueue.emplace_back(l,n,a,lprev,nprev,aprev);
        }         
        std::cout << "FFSConfigIDQueue has " << FFSConfigIDQueue.size() << " entries upon initialization\n";

        // now that queue is populated initialize tasks for all processors
        // ==============================

        bool shouldpop_local = true
        PopQueueMPI(sholdpop_local);
    }

    void ForwardFlux::PrintQueue(){
        for (int i =0 ;i < FFSConfigIDQueue.size(); i++){
            std::cout << i <<" "
                      <<FFSConfigIDQueue[i].l <<" "
                      <<FFSConfigIDQueue[i].n <<" "
                      <<FFSConfigIDQueue[i].a <<" "
                      <<FFSConfigIDQueue[i].lprev <<" "
                      <<FFSConfigIDQueue[i].nprev <<" "
                      <<FFSConfigIDQueue[i].aprev <<"\n";
        }
    }

    void ForwardFlux::PopQueueMPI(bool shouldpop_local){

        bool *shouldpop = new bool(_world.size());

        MPI_Allgather(&shouldpop_local,1,MPI::BOOL,shouldpop,1,MPI::BOOL,_world);

        // I dont pass the queue information between procs but I do syncronize 'shouldpop'
        //   as a reuslt all proc should have the same queue throughout the simulation
        for (int i=0;i<_world.size();i++){
          if (shouldpop[i] == true){ 
            if (i == _world.rank()){ //if rank matches read and pop
              if (!FFSConfigIDQueue.empty()){ //if queue has tasks
                   myFFSConfigID = FFSConfigIDQueue.front();
                   ReadFFSConfiguration (snapshot, myFFSConfigID);
                   _pop_tried_but_empty_queue = false;
                  FFSConfigIDQueue.pop_front();
              }
              else{ //queue is empty, need to wait for new tasks to come in
                if ((successes[i] == true) || (failures[i] == true)) _pop_tried_but_empty_queue = true;
              }
            }
            else{ //else if rank doesnt match, just pop 
              if (!FFSConfigIDQueue.empty()){
                  FFSConfigIDQueue.pop_front();
              }
            }
          }
        }
        // done copy (perhaps I should make this a method)
        // ==============================



    }
    

	
}

