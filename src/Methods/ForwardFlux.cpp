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
        _ninterfaces = 5;
        int i;
        _interfaces.resize(_ninterfaces);
        _interfaces[0]=-1.0;
        _interfaces[1]=-0.8;
        _interfaces[2]=-0.6;
        _interfaces[3]= 0;
        _interfaces[4]= 1.0;

        _saveTrajectories = true;

        _current_interface = 0;

        _output_directory = "FFSoutput";
        _initialFluxFlag = true;
        _M.resize(_ninterfaces);
        _A.resize(_ninterfaces);
        _P.resize(_ninterfaces);
        _S.resize(_ninterfaces);
        _N.resize(_ninterfaces);
        for(i=0;i<_ninterfaces;i++) _M[i] = 50;

        _N[0] = 100;
        
        Lambda0ConfigLibrary.resize(_N[0]);
        std::normal_distribution<double> distribution(0,1);
        for (i = 0; i < _N[0] ; i++){
          Lambda0ConfigLibrary[i].l = 0;
          Lambda0ConfigLibrary[i].n = i;
          Lambda0ConfigLibrary[i].a = 0;
          Lambda0ConfigLibrary[i].lprev = 0;
          Lambda0ConfigLibrary[i].nprev = i;
          Lambda0ConfigLibrary[i].aprev = 0;

        // I beleive the PreSimulation must run once at the beginning of the simulation. Apparently it runs every iteration!!!
        
          /*FFSConfigID ffsconfig = Lambda0ConfigLibrary[i];
          // Write the dump file out
          std::ofstream file;
          std::string filename = _output_directory + "/l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + ".dat";
          file.open(filename.c_str());

          //first line gives ID of where it came from
          file << ffsconfig.lprev << " " << ffsconfig.nprev << " " << ffsconfig.aprev << "\n";
          //write position and velocity
          file << "1 -1 0 0 " << distribution(_generator) << " "  << distribution(_generator) << " 0\n";*/

        }
        _pop_tried_but_empty_queue = false;


    }

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
        //check if we want to check FFS interfaces this timestep
        //for now, do it every time step
        if (_iteration % 1 != 0) return;

        // check the structure at the beginning of the simulation
        if (_iteration == 0) {
          CheckInitialStructure(cvs);
        }

        // if _computefluxA0
        if (_initialFluxFlag){
            ComputeInitialFlux(snapshot,cvs); 
            if (!_initialFluxFlag){ //only enter here once

              InitializeQueue(snapshot,cvs);
              PrintQueue();
            }
        }
        // Else check the FFS interfaces
        else{
            CheckForInterfaceCrossings(snapshot,cvs);
            //FluxBruteForce(snapshot,cvs);

        }
        // Other modes?

    }

	void ForwardFlux::PostSimulation(Snapshot* snapshot, const CVList& cvs){

        if (_saveTrajectories){
          ReconstructTrajectories(snapshot);	
        }

        _world.abort(EXIT_FAILURE); //more elegant solution?
	}

  void ForwardFlux::CheckInitialStructure(const CVList& cvs){

        if (_initialFluxFlag){
          std::cout << "Running initial Flux calculations" << std::endl;

          // Check if we are in State A
          _cvvalue = cvs[0]->GetValue();
          double _firstInterfaceLocation = _interfaces[0];
          if ( _cvvalue > _firstInterfaceLocation) {
            std::cout << "Please provide an initial configuration in State A. Exiting ...." << std::endl;
            exit(1);
          }
        }
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

        unsigned int nInitialConfigs = _N[0];

        bool shouldContinueLocal = true;
        InitialFluxMPI(snapshot, cvs, shouldContinueLocal);

        //check if we've crossed the first interface (lambda 0)
        int hascrossed = HasCrossedInterface(_cvvalue, _cvvalue_previous, 0);
        bool success_local = false;
        bool *successes = new bool (_world.size());

        // If we have crossed the first interface, then write the information to disk
        if (hascrossed == 1){
              success_local = true;
              //for each traj that crossed to lambda0 in forward direction, we need to write it to disk (FFSConfigurationFile)
              MPI_Allgather(&success_local,1,MPI::BOOL,successes,1,MPI::BOOL,_world);

              for (int i = 0; i < _world.size(); i++){
                int l,n,a,lprev,nprev,aprev;
                // Since we are in State A, the values of lprev, nprev, aprev are all zero.
                lprev = 0;
                nprev = 0;
                aprev = 0;

                if (successes[i] == true){ 
                  if (i == _world.rank()){
                    //update ffsconfigid's l,n,a
                    l = 0;
                    n = _S[0];
                    a = 0;
                    FFSConfigID newid = FFSConfigID(l,n,a,lprev,nprev,aprev);
                    WriteFFSConfiguration(snapshot,newid,1);
                    _S[0]++;
                  }
                }    
              }  
            }

        //print some info
        if (shouldContinueLocal){
          std::cout << "Iteration: "<< _iteration << ", proc " << _world.rank() << std::endl;
          if (success_local)
            std::cout << "Successful attempt. (cvvalue_previous: " << _cvvalue_previous << " cvvalue " << _cvvalue << " )" << std::endl;
          std::cout << "# of successes:               " << _S[0] << std::endl;
          std::cout << "required # of configurations: " << nInitialConfigs << std::endl;
        }

        // Check if the required number of initial configurations are created, if so print a message, and compute the initial flux.
        if (_S[0] == nInitialConfigs){
            std::cout << "Initial flux calculation was successfully completed" << std::endl;
            shouldContinueLocal = false;            
            // Call a function to compute the actual flux and output the data
            WriteInitialFlux(snapshot, cvs);
            // set _S[0] to zero
            _S[0] = 0;
            //responsible for setting _initialFluxFlag = false when finished
            _initialFluxFlag = false;
            //exit(1);
          } 

        _cvvalue_previous = _cvvalue;
        _iteration++;

        //clean up
        delete[] successes;
    }

  void ForwardFlux::InitialFluxMPI(Snapshot* snapshot, const CVList& cvs, bool shouldContinueLocal){

      bool *shouldContinue = new bool(_world.size());
      MPI_Allgather(&shouldContinueLocal, 1, MPI::BOOL, shouldContinue, 1, MPI::BOOL, _world);

      for (int i = 0; i < _world.size(); i++){
        if (shouldContinue[i] == true){
          if(i == _world.rank()) {
            // Trigger a rebuild of the cvs since we reset the positions
              cvs[0]->Evaluate(*snapshot);
              _cvvalue = cvs[0]->GetValue();
          }
        }
      }

      //clean up
      delete[] shouldContinue;
    }

  void ForwardFlux::WriteInitialFlux(Snapshot* snapshot, const CVList& cvs){

      std::ofstream file;
      std::string filename = _output_directory + "/initial_flux_value.dat";
      file.open(filename.c_str());
      if (!file){
        std::cerr << "Error! Unable to write " << filename << std::endl;
        exit(1);
      }
      _N0TotalSimTime = _iteration * _world.size();
      _fluxA0 = (double) (_N[0] / _N0TotalSimTime);
      file << "number of processors: " << _world.size() << std::endl;
      file << "number of iterations: " << _iteration << std::endl;
      file << "Total simulation time: " << _N0TotalSimTime << std::endl;
      file << "Initial flux: " << _fluxA0 << std::endl;
      file.close();
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
        //std::vector<bool> in MPI were strange, ended up using arrays
        bool *successes = new bool (_world.size());
        bool *failures = new bool (_world.size());

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
          int l,n,a,lprev,nprev,aprev;
          // write config to lambda+1
          lprev = myFFSConfigID.l;
          nprev = myFFSConfigID.n;
          aprev = myFFSConfigID.a;

          if (successes[i] == true){ 
            if (i == _world.rank()){
              //update ffsconfigid's l,n,a
              l = _current_interface + 1;
              n = _S[_current_interface] + success_count;
              a = 0;
              FFSConfigID newid = FFSConfigID(l,n,a,lprev,nprev,aprev);
              WriteFFSConfiguration(snapshot,newid,1);
            }
            success_count++;
          }
          if (failures[i] == true){ 
            if (i == _world.rank()){
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
        
        //update the number of successes and attempts, same for all proc since Allgathered 'successes' and 'failures'
        _S[_current_interface] += success_count;
        _A[_current_interface] += success_count + fail_count;
        // ^ I dont like storing attempts this way (as is its only when they finish). _A should also include jobs in progress (i.e. jobs currently running). THINK ABOUT THIS!.
        _nfailure_total += fail_count;
        
        // Check if this interface is finished, if so add new tasks to queue, and increment _current_interface
        if (_A[_current_interface] == _M[_current_interface]-1){
          if (_current_interface+1 != _ninterfaces){
            _current_interface += 1;
            _N[_current_interface] = _S[_current_interface-1];

            if (_N[_current_interface] == 0){
               std::cerr << "Error! No successes from interface " << _current_interface-1 << " to " << _current_interface <<"! Try turning up M["<<_current_interface-1<<"] or spacing interfaces closer together.\n";
              _world.abort(EXIT_FAILURE); 
            }

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
              l = _current_interface;
              n = mypick;
              a = attempt_count[mypick]; 
              attempt_count[mypick]++; //this updates attempt number if same config is picked twice

              FFSConfigIDQueue.emplace_back(l,n,a,lprev,nprev,aprev);
            }
          }
          else{
            std::cout << "DFFS should be finished here, do something special? like exit?\n";
            //Hythem said this is "acceptable" until the code is changed
            PostSimulation(snapshot,cvs);
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
        
        if (_saveTrajectories){
          if (!_pop_tried_but_empty_queue){ //dont update trajectory if zombie job
              AppendTrajectoryFile(snapshot,_trajectory_file);
          }
          if (success_local || fail_local){ //if finished then close it
            _trajectory_file.close();
          }
        }
        //Pop the queue
        // Need to perform mpi call so that all proc pop the queue in the same way
        PopQueueMPI(snapshot,cvs, shouldpop_local);
                
        //Anything else to update across mpi?

        //print info
        if (shouldpop_local){
          std::cout << "Iteration: "<< _iteration << ", proc " << _world.rank() << "\n";
          if (success_local)
            std::cout << "Successful attempt from interface " << _current_interface <<" (cvvalue_previous: " << _cvvalue_previous << " cvvalue " << _cvvalue << " interface " << _interfaces[_current_interface] << "\n";
          if (fail_local)
            std::cout << "Failed attempt from interface " << _current_interface <<" (cvvalue_previous: " << _cvvalue_previous << " cvvalue " << _cvvalue << " interface " << _interfaces[0] << "\n";
          std::cout << "A: ";
          for (auto a : _A) std::cout << a << " "; std::cout << "\n";
          std::cout << "S: ";
          for (auto s : _S) std::cout << s << " "; std::cout << "\n";
          std::cout << "M: ";
          for (auto m : _M) std::cout << m << " "; std::cout << "\n";
        }
        
        //clean up
        _cvvalue_previous = _cvvalue;
        _iteration++;

        delete[] successes;
        delete[] failures,shouldpop;
    }

	void ForwardFlux::ComputeTransitionProbabilities(){

    }

	
	void ForwardFlux::WriteFFSConfiguration(Snapshot* snapshot, FFSConfigID& ffsconfig, bool wassuccess)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();
		//const auto& dumpfilename = snapshot->GetSnapshotID();

        // Write the dump file out
		std::ofstream file;
        std::string filename;
        if (wassuccess){
          filename = _output_directory + "/l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + ".dat";
        }
        else{
          filename = _output_directory + "/fail-l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + ".dat";

        }
 		file.open(filename.c_str());
        if (!file) {std::cerr << "Error! Unable to write " << filename << "\n"; exit(1);}

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

    void ForwardFlux::OpenTrajectoryFile(std::ofstream& file){

        FFSConfigID ffsconfig = myFFSConfigID;

        std::string filename = _output_directory + "/traj-l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + "-a" + std::to_string(ffsconfig.a) + ".xyz";
 		file.open(filename.c_str());
        if (!file) {std::cerr << "Error! Unable to write " << filename << "\n"; exit(1);}
    
    }

    void ForwardFlux::AppendTrajectoryFile(Snapshot* snapshot, std::ofstream& file){

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

	void ForwardFlux::ReadFFSConfiguration(Snapshot* snapshot, FFSConfigID& ffsconfig)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		//auto& ID = snapshot->GetSnapshotID();

        std::ifstream file; 
        std::string filename =  _output_directory + "/l"+ std::to_string(ffsconfig.l) +"-n"+ std::to_string(ffsconfig.n) + ".dat";
        std::string line;

        file.open(filename);
        if (!file) {std::cerr << "Error! Unable to read " << filename << "\n"; exit(1);}

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
          std::uniform_int_distribution<int> distribution(0,_N[0]-1);
          for (int i=0; i < npicks ; i++){
             picks[i] = distribution(_generator);
          }
        }
        MPI_Bcast(picks.data(),npicks,MPI::UNSIGNED,0,_world);

        //set correct attempt index if a given ID is picked twice
        std::vector<unsigned int> attempt_count;
        attempt_count.resize(_N[0],0);

        //each proc adds to the queue
        for (int i=0; i < npicks ; i++){
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

        // now that queue is populated initialize tasks for all processors
        // ==============================

        bool shouldpop_local = true;
        
        PopQueueMPI(snapshot,cvs,shouldpop_local);

        
        _cvvalue_previous = cvs[0]->GetValue();

    }

    void ForwardFlux::PrintQueue(){
        for (unsigned int i =0 ;i < FFSConfigIDQueue.size(); i++){
            std::cout << i <<" "
                      <<FFSConfigIDQueue[i].l <<" "
                      <<FFSConfigIDQueue[i].n <<" "
                      <<FFSConfigIDQueue[i].a <<" "
                      <<FFSConfigIDQueue[i].lprev <<" "
                      <<FFSConfigIDQueue[i].nprev <<" "
                      <<FFSConfigIDQueue[i].aprev <<"\n";
        }
    }

    void ForwardFlux::PopQueueMPI(Snapshot* snapshot, const CVList& cvs, bool shouldpop_local){

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
                   
                   //open new trajectory file, write first frame
                   if (_saveTrajectories){ 
                     OpenTrajectoryFile(_trajectory_file);
                     AppendTrajectoryFile(snapshot,_trajectory_file);
                   }

                   //Trigger a rebuild of the CVs since we reset the positions
                   cvs[0]->Evaluate(*snapshot);
                   _cvvalue_previous = cvs[0]->GetValue();

                  _pop_tried_but_empty_queue = false;
                  FFSConfigIDQueue.pop_front();
              }
              else{ //queue is empty, need to wait for new tasks to come in
                 _pop_tried_but_empty_queue = true;
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
    
    void ForwardFlux::ReconstructTrajectories(Snapshot *snapshot){
        int nsuccess = _S[_ninterfaces-1];
        
        
        //Snapshot* snapshot_empty;
        FFSConfigID ffsconfig;

        for(int i = 0; i < nsuccess ; i++){
          int lprev,nprev,aprev;
          std::deque<FFSConfigID> path;

          ffsconfig.l = _ninterfaces-1;
          ffsconfig.n = i;
          ffsconfig.a = 0;
          while (ffsconfig.l != 0){
             //this is just to populate ffsconfig.{lprev,nprev,aprev}
             //ReadFFSConfiguration(snapshot_empty,ffsconfig);
             ReadFFSConfiguration(snapshot,ffsconfig);
             
             //path.emplace_front(ffsconfig); //not sure if new constructor will work
             path.emplace_front(ffsconfig.l,ffsconfig.n, ffsconfig.a, ffsconfig.lprev, ffsconfig.nprev, ffsconfig.aprev); 

             ffsconfig.l = ffsconfig.lprev;
             ffsconfig.n = ffsconfig.nprev;
             ffsconfig.a = ffsconfig.aprev;
          }

          //now path should contain all of the FFSConfigID's from B back to A
          //reverse pop it and splice all traj- files into a new traj-full- file

          //output file
          std::ofstream ofile;
          std::string ofilename;
          ofilename = _output_directory +"/traj-full-" + std::to_string(i) + ".xyz";
          ofile.open(ofilename.c_str());
          if (!ofile) {std::cerr << "Error! Unable to write " << ofilename << "\n"; exit(1);}

          while(!path.empty()){
            ffsconfig = path.front();
            path.pop_front();

            //input file
            std::ifstream ifile;
            std::string ifilename;

            ifilename = _output_directory + "/traj-l" + std::to_string(ffsconfig.l) + "-n" + std::to_string(ffsconfig.n) + "-a" + std::to_string(ffsconfig.a) + ".xyz";

            ifile.open(ofilename.c_str());
            if (!ifile) {std::cerr << "Error! Unable to read " << ifilename << "\n"; exit(1);}
            //write entire ifile to ofile
            std::string line;
            while(!std::getline(ifile,line).eof()){
                ofile << line;
            }
          }
          ofile.close();
        }
    }
    

    void ForwardFlux::FluxBruteForce(Snapshot* snapshot, const CVList& cvs){
        //run a simulation, if it returns to A get a new config, if it suceeds, store the time it took

        //eventually will overwrite CheckForInterfaceCrossings (or copy from it)
        // I want to make sure everything is correct with FFS first though
        


    }
	
}

