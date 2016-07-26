#include "Swarm.h"
#include "../spline.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>

namespace mpi = boost::mpi;
namespace SSAGES
{

    //Helper function for calculating distances
    double Swarm::distance(std::vector<double>& x, std::vector<double>& y)
    {
        double distance = 0;
        for(size_t i = 0; i < x.size(); i++)
        {
            //Calculate the distance between each CV in x and y
            distance += pow((x[i]-y[i]),2.0);
        }
        return sqrtl(distance);
    }

    //Helper function to check if CVs are initialized correctly
    bool Swarm::CVInitialized(const CVList& cvs)
    {
        //std::cout << _mpiid << " Made a call to CVInitialized()" << std::endl; //Debugging
        double threshold = 0.05;
        const double epsilon = 0.000000001;
        double diff;

        //On the first iteration, check that the CVs are within (threshold*100)% of the center value they're associated to
        for(size_t i = 0; i < cvs.size(); i++)
        { 
            std::cout << _mpiid << " CV value = " << cvs[i]->GetValue() << " Centers value = " << _centers[i] << std::endl; //Debuggingi
            if(_centers[i] <= epsilon)
            {//e.g. if _centers[i] = 0
                diff = std::abs((cvs[i]->GetValue() - (_centers[i]+0.01)) / ((cvs[i]->GetValue() + (0.01 + _centers[i]))/2.0));
            }
            else
            {
                diff = std::abs((cvs[i]->GetValue() - (_centers[i])) / ((cvs[i]->GetValue() + (_centers[i]))/2.0));
            }
            //std::cout << _mpiid << " Diff = " << diff << "Abs Test = " << std::abs(-10.1) << std::endl; //Debugging
            if(diff >= threshold)
            {
                std::cout << _mpiid << " Diff true" << std::endl; //Debugging
                return true; //e.g. proceed to initialize again
            }
        }
        //std::cout << _mpiid << " Diff false" << std::endl; //Debugging
        return false; //e.g. OK to move on to regular sampling
    }

    //Pre-simulation hook
    void Swarm::PreSimulation(Snapshot* snapshot, const CVList& cvs)
    {
        auto& forces = snapshot->GetForces();
        auto& positions = snapshot->GetPositions();
        auto& velocities = snapshot->GetVelocities();

        //Open file for writing
        _mpiid = snapshot->GetWalkerID();
        std::ostringstream file;
        file << "node-" << std::setw(4) << std::setfill('0') << _mpiid << ".log";
        _stringout.open(file.str());
        //Sizing vectors
        _worldstring.resize(_centers.size()); 
        _cv_start.resize(_centers.size());
        _cv_drift.resize(_centers.size());
        _traj_positions.resize(_number_trajectories);
        _traj_forces.resize(_number_trajectories);
        _traj_velocities.resize(_number_trajectories);

        //std::cout << _mpiid <<  " Reserving size..." << std::endl; //Debugging
        _traj_positions.reserve(_number_trajectories * positions.size());
        _traj_forces.reserve(_number_trajectories * forces.size());
        _traj_velocities.reserve(_number_trajectories * velocities.size());
        
        for(size_t k = 0; k < _traj_positions.size(); k++)
        {
            _traj_positions[k].resize(positions.size());
            _traj_positions[k].reserve(positions.size());
        }
        for(size_t k = 0; k < _traj_forces.size(); k++)
        {
            _traj_forces[k].resize(forces.size());
            _traj_forces[k].reserve(forces.size());
        }
        for(size_t k = 0; k < _traj_velocities.size(); k++)
        {
            _traj_velocities[k].resize(velocities.size());
            _traj_velocities[k].reserve(velocities.size());
        }

        //Initialize vector values
        for(size_t i = 0; i < _centers.size(); i++)
        {
            _worldstring[i].resize(_numnodes); 

            //_worldstring is indexed as cv followed by node
            mpi::all_gather(_world, _centers[i], _worldstring[i]);

            _cv_start[i] = 0;
            _cv_drift[i] = 0; 
        }

        //Additional initializing
        _index = 0;  
        _restrained_steps = _harvest_length*_number_trajectories; 
        _unrestrained_steps = _swarm_length*_number_trajectories;

        //For reparametrization
        _alpha = _mpiid / (_numnodes - 1.0);

        PrintString(cvs);

        //Print initial string
        if(_mpiid == 0)
        {
            std::ofstream pyth_string ("2dstring_init.txt");
            for(size_t i = 0; i < _numnodes; i++)
            {
                for(size_t j = 0; j < _centers.size(); j++)
                {
                    pyth_string << _worldstring[j][i] << " ";
                    if(j == _centers.size()-1)
                    {
                        pyth_string << std::endl;
                    }
                }
            }
            pyth_string.close();
        }

        sampling_started = true;
        IntegrationCounter = 0;
    }

    void Swarm::PostIntegration(Snapshot* snapshot, const CVList& cvs)
    {
        auto& forces = snapshot->GetForces();
        auto& positions = snapshot->GetPositions();
        auto& velocities = snapshot->GetVelocities();

        bool initialize; //Whether to initialize or not

        std::cout << _mpiid << " PostIntegration, Number " << IntegrationCounter << std::endl;
        IntegrationCounter++;

        if(_currentiter == 0 && sampling_started) 
        {
            initialize = CVInitialized(cvs);
        }
        else
        {
            initialize = false;
        }
        if(_currentiter == 0 && initialize && sampling_started)
        {//On first pass, make sure CVs are initialized well
            //std::cout << _mpiid << " Initializing on first pass..." << std::endl; //Debugging
            //Do restrained sampling, and do not harvest trajectories
            for(size_t i = 0; i < cvs.size(); i++)
            {
                std::cout << _mpiid << " First Pass Restraining..." << std::endl; //Debugging
                //Get current CV and gradient
                auto& cv = cvs[i];
                auto& grad = cv->GetGradient();

                //Compute dV/dCV
                auto D = _spring*(cv->GetDifference(_centers[i]));

                //Update forces
                for(size_t j = 0; j < forces.size(); j++)
                {
                    for(size_t k = 0; k < forces[j].size(); k++)
                    {
                        forces[j][k] -= (double)D*grad[j][k]; 
                        //forces[j][k] -= 0.0; //Debugging
                    }
                }
            }
        }
        else
        {
            if(sampling_started)
            {
                sampling_started = false; //Flag to prevent unneeded umbrella sampling
            }
            if(_iterator <= _initialize_steps + _restrained_steps)
            {
                //Do restrained sampling, and do not harvest trajectories
                for(size_t i = 0; i < cvs.size(); i++)
                {
                    std::cout << _mpiid << " Normal  Restraining..." << std::endl; //Debugging
                    if(_iterator == 0)
                    {
                        _index = 0; //Reset index when starting
                    }
                    //Get current CV and gradient
                    auto& cv = cvs[i];
                    auto& grad = cv->GetGradient();

                    //Compute dV/dCV
                    auto D = _spring*(cv->GetDifference(_centers[i]));

                    //Update forces
                    for(size_t j = 0; j < forces.size(); j++)
                    {
                        for(size_t k = 0; k < forces[j].size(); k++)
                        {
                            forces[j][k] -= (double)D*grad[j][k]; 
                            //forces[j][k] -= 0.0; //Debugging

                        }
                    }
                }
                if(_iterator > _initialize_steps)
                {
                    //Harvest a trajectory every ten steps
                    if(_iterator % 10 == 0)
                    {
                        std::cout << _mpiid << " Harvesting" << std::endl; //Debugging
                        for(size_t k = 0; k < positions.size(); k++)
                        {
                            for(size_t l = 0; l < positions[k].size(); l++)
                            {
                                _traj_positions[_index][k][l] = positions[k][l];
                            }
                        }
                        for(size_t k = 0; k < forces.size(); k++)
                        {
                            for(size_t l = 0; l < forces[k].size(); l++)
                            {
                                _traj_forces[_index][k][l] = forces[k][l];
                            }
                        }
                        for(size_t k = 0; k < velocities.size(); k++)
                        {
                            for(size_t l = 0; l < velocities[k].size(); l++)
                            {
                                _traj_velocities[_index][k][l] = velocities[k][l];
                            }
                        }
                        _index++;
                    }
                }
                if(_iterator == _initialize_steps + _restrained_steps)
                {
                    //Reset positions and forces before first call to unrestrained sampling
                    _index = 0;
                    {
                        std::cout << _mpiid << " Attempting reset" << std::endl;
                        //Zero forces
                        for(auto& force: forces)
                            for(auto& xyz : force)
                                xyz = 0.0;
                        //Then set positions
                        for(size_t k = 0; k < positions.size(); k++)
                        {
                            for(size_t l = 0; l < positions[k].size(); l++)
                            {
                                positions[k][l] = _traj_positions[_index][k][l];
                            }
                        }
                        //Then set forces
                        for(size_t k = 0; k < forces.size(); k++)
                        {
                            for(size_t l = 0; l < forces[k].size(); l++)
                            {
                                //forces[k][l] = _traj_forces[_index][k][l];
                                forces[k][l] = 0.0;
                            }
                        }
                        //Velocities
                        for(size_t k = 0; k < velocities.size(); k++)
                        {
                            for(size_t l = 0; l < velocities[k].size(); l++)
                            {
                                velocities[k][l] = _traj_velocities[_index][k][l];
                            }
                        }
                    }
                }
                _iterator++;
            }
            else if(_iterator <= _initialize_steps + _restrained_steps + _unrestrained_steps)
            {
                //Launch unrestrained trajectories
                //std::cout << _mpiid <<  " Running swarm...Iteration number = " << _iterator << std::endl; //Debugging  
                if((_iterator - _initialize_steps - _restrained_steps) % _swarm_length == 0)
                {
                    //std::cout << _mpiid << " End of trajectory..." << std::endl; //Debugging
                    //End of trajectory, harvest drift
                    for(size_t i = 0; i < _cv_drift.size(); i++)
                    {
                        //Drift scaled down for debugging
                        _cv_drift[i] = (_cv_drift[i]*_index + cvs[i]->GetValue()  - _centers[i]) / (_index+1); //Calculate running average of drifts
                        //_cv_drift[i] += 0; //Debugging
                        //std::cout << _mpiid << " CV Drift = " << _cv_drift[i] << " " << std::endl; //Debugging
                    }
                    //Set up for next trajectory
                    _index++;
                    if(_index < _number_trajectories)
                    {
                        //std::cout << _mpiid << " Starting trajectory...Index == " << _index << std::endl; //Debugging
                        //std::cout << _mpiid << " Current size of forces and positions " << _traj_forces.size() << " " << _traj_positions.size() << std::endl;
                        //Start of trajectory, reset positions and forces
                        std::cout << _mpiid << " Attempting reset" << std::endl;
                        {
                            //Debugging?  Start be zeroing
                            for(auto& force : forces)
                                for(auto& xyz : force)
                                    xyz = 0.0;
                            //Then set positions
                            for(size_t k = 0; k < positions.size(); k++)
                            {
                                for(size_t l = 0; l < positions[k].size(); l++)
                                {
                                    positions[k][l] = _traj_positions[_index][k][l];
                                }
                            }
                            //Then set forces
                            for(size_t k = 0; k < forces.size(); k++)
                            {
                                for(size_t l = 0; l < forces[k].size(); l++)
                                {
                                    forces[k][l] = 0.0;
                                    //forces[k][l] = _traj_forces[_index][k][l];
                                }
                            }
                            //Velocities
                            for(size_t k = 0; k < velocities.size(); k++)
                            {
                                for(size_t l = 0; l < velocities[k].size(); l++)
                                {
                                    velocities[k][l] = _traj_velocities[_index][k][l];
                                }
                            }   
                        }
                    }
                }
                _iterator++;

                if(_iterator == _initialize_steps + _restrained_steps + _unrestrained_steps + 1)
                {
                    //std::cout << _mpiid << " Last trajectory call" << std::endl; //Debugging
                }
            }
            else
            {
                //std::cout << _mpiid << " Accessed final loop" << std::endl; //Debugging
                //_world.barrier(); //Hold until everything gets here
                //std::cout << _mpiid << " Starting CV update" << std::endl; //Debugging
                //Average drift
                /*for(size_t i = 0; i < _cv_drift.size(); i++)
                {
                    _cv_drift[i] /= _number_trajectories;
                }*/

                //Evolve CVs, reparametrize, and reset vectors
                _currentiter++;
                //std::cout << _mpiid << " Reached string iteration " << _currentiter << std::endl; //Debugging

                _world.barrier(); //Wait for all nodes before attempting string update
                StringUpdate();
                _world.barrier(); //Wait for all CVs to be updated
                PrintString(cvs);

                _iterator = 0;
                _index = 0;

                for(size_t i = 0; i < _cv_drift.size(); i++)
                {
                    _cv_drift[i] = 0; 
                }
                //_world.barrier(); //Hold until all CVs are updated
            }
        }
    }

    //Post simulation hook
    void Swarm::PostSimulation(Snapshot*, const CVList&)
    {
        _stringout.close();
    }

    void Swarm::PrintString(const CVList& CV)
    {
        std::cout << _mpiid << " Printing string" << std::endl; //Debugging
        //Write node, iteration, centers of the string and current CV value to output file
        _stringout.precision(8);
        _stringout << _mpiid << " " << _currentiter << " ";

        for(size_t i = 0; i < _centers.size(); i++)
        {
            _stringout << _centers[i] << " " << CV[i]->GetValue() << " ";
        }
        _stringout << std::endl;

        //Write same info to terminal, omit current CV value
        std::cout << _mpiid << " " << _currentiter << " "; 
        for(size_t i = 0; i < _centers.size(); i++)
        {
            std::cout << _centers[i] << " "; 
        }
        std::cout << std::endl;
    }

    void Swarm::StringUpdate()
    {
        std::cout << _mpiid << " Updating string..." << std::endl; //Debugging
        //Values for reparametrization
        size_t i;
        int centersize = _centers.size();
        double alpha_star;
        std::vector<double> alpha_star_vector;
        std::vector<double> cvs_new; 
        std::vector<double> cvs_new_vector;

        //MPI communication parameters
        int sendneighbor, recvneighbor;
        MPI_Status status;

        //Vectors for lower and upper neighbor on the string
        std::vector<double> lower_cv_neighbor;
        lower_cv_neighbor.resize(centersize); 

        //For reparametrization
        cvs_new.resize(centersize); //All CVs for a particular nodes
        alpha_star_vector.resize(_numnodes);
        cvs_new_vector.resize(_numnodes); //Eventually, all CVs of one dimension from each node

        //Evolve CVs with average drift
        for(i = 0; i < cvs_new.size(); i++)
        {
            cvs_new[i] = _centers[i] + _cv_drift[i]*_drift_scale; //Rescale drift
            //std::cout << cvs_new[i] << " " << std::endl; //Debugging
        }

        //Set up nodes to receive from their backward neighbor and send to their forward neighbor, wrapping around at the string ends
        if(_mpiid == 0)
        {
            sendneighbor = 1;
            recvneighbor = _world.size() - 1;
        }
        else if (_mpiid == _world.size() - 1)
        {
            sendneighbor = 0;
            recvneighbor = _world.size() - 2;
        }
        else
        {
            sendneighbor = _mpiid + 1;
            recvneighbor = _mpiid - 1;
        }

        //Copy centers from lower neighbor into lower_cv_neighbor
        MPI_Sendrecv(&_centers[0], centersize, MPI_DOUBLE, sendneighbor, 1234, &lower_cv_neighbor[0], centersize, MPI_DOUBLE, recvneighbor, 1234, _world, &status);

        //Reparameterization
        //Alpha star is the uneven mesh, approximated as linear distance between points
        if(_mpiid == 0)
        {
            alpha_star = 0;
        }
        else
        {
            alpha_star = distance(_centers, lower_cv_neighbor);
        }

        //Gather each alpha_star into a vector 
        mpi::all_gather(_world, alpha_star, alpha_star_vector);

        for(i = 1; i < alpha_star_vector.size(); i++)
        {
            //Distances should be summed
            alpha_star_vector[i] = alpha_star_vector[i-1] + alpha_star_vector[i];
        }
        
        for(i = 1; i < alpha_star_vector.size(); i++)
        {
            //Distances should be normalized
            alpha_star_vector[i] /= alpha_star_vector[_numnodes - 1];
        }

        tk::spline spl; //Cubic spline interpolation

        for(i = 0; i < centersize; i++)
        {
            //Take from new CV values in one dimension into a vector
            //cvs_new_vector is thus the CV value in one dimension at each node
            mpi::all_gather(_world, cvs_new[i], cvs_new_vector);
            //spl.setpoints(alpha_star_vector, vector of all new points in one particular dimension)
            spl.set_points(alpha_star_vector, cvs_new_vector); //Spline initialized with points on uneven mesh
            //Update CV values to lie on a regular mesh
            _centers[i] = spl(_alpha); 
        }
        //Print a simple text file to be analyzed by Python script
        for(size_t i = 0; i < _centers.size(); i++)
        {//Fill up world_string
            //_worldstring is indexed as cv followed by node
            mpi::all_gather(_world, _centers[i], _worldstring[i]);
        }
        //_world.barrier();
        if(_mpiid == 0)
        {
            std::ofstream pyth_string ("2dstring.txt");
            for(size_t i = 0; i < _numnodes; i++)
            {
                for(size_t j = 0; j < centersize; j++)
                {
                    pyth_string << _worldstring[j][i] << " ";
                    if(j == centersize-1)
                    {
                        pyth_string << std::endl;
                    }
                }
            }
            pyth_string.close();
        }
    }
}

