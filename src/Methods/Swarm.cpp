#include "Swarm.h"
#include "../spline.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace mpi = boost::mpi;
namespace SSAGES
{

    //Helper function for calculating distances
    double distance(std::vector<double>& x, std::vector<double>& y)
    {
        double distance = 0;
        for(size_t i = 0; i < x.size(); i++)
        {
            //Calculate the distance between each CV in x and y
            distance += pow((x[i]-y[i]),2.0);
        }
        return sqrtl(distance);
    }

    //Pre-simulation hook
    void Swarm::PreSimuation(Snapshot* snapshot, const CVList& cvs)
    {
        //Open file for writing
        _mpiid = snapshot->GetWalkerID();
        auto& positions = snapshot->GetPositions();
        char file[1024];
        sprintf(file, "node-%04.d.log", _mpiid);
        _stringout.open(file);

        //Sizing vectors
        _worldstring.resize(_centers.size()); 
        _cv_start.resize(_centers.size();
        _cv_drift.resize(_centers.size());
        _traj_positions.resize(_number_trajectories);
        _traj_forces.resize(_number_trajectories);

        //Initialize vector values
        for(size_t i = 0; i < _centers.size(); i++)
        {
            _worldstring[i].resize(_numnodes);
            _cv_prev[i] = cvs[i]->GetValue();

            //_worldstring is indexed as cv followed by node
            mpi::all_gather(_world, _centers[i], _worldstring[i]);

            _cv_start[i] = 0;
            _cv_drift[i] = 0; 
        }

        //Additional initializing
        index = 0;  
        _restrained_steps = _harvest_length*_number_trajectories; 
        _unrestrained_steps = _swarm_length*number_trajectories;
    }

    void Swarm::PostIntegration(Snapshot* snapshot, const CVList& cvs)
    {
        auto& forces = snapshot->GetForces();
        auto& positions = snapshot->GetPositions();
        int index; //For building the trajectory vectors

        if(_iterator <= _initialize_steps + _restrained_steps)
        {
            //Do restrained sampling, and do not harvest trajectories
            for(size_t i = 0; i < cvs.size(); i++)
            {
                if(_iterator = 0)
                {
                    index = 0; //Reset index when starting
                }
                //Get current CV and gradient
                auto& cv = cvs[i];
                auto& grad  cv->GetGradient();

                //Compute dV/dCV
                auto D = _spring*(cv->GetDifference(_centers[i]));

                //Update forces
                for(size_t j = 0; j < forces.size(), j++)
                {
                    for(size_t k = 0; k < forces[j].size(), k++)
                    {
                        forces[j][k] -= D*grad[j][k]; 

                    }
                }
            }
            if(_iterator > _initialize_steps)
            {
                //Harvest a trajectory every ten steps
                if(_iterator % 10 == 0)
                {
                    _traj_positions[index] = positions;
                    _traj_forces[index] = forces; 
                    index++;
                }
            }
            _iterator++;
        }
        else if(_iterator <= _initialize_steps + _restrained_steps + _unrestrained_steps)
        {
            //Launch unrestrained trajectories
            if(_iterator == _initialize_steps + _restrained_steps + 1)
            {
                index = 0; //Reset index when starting unrestrained trajectories
            }
            if((iterator - _initialize_steps - restrained_steps) % swarm_length == 1)
            {
                //Start of trajectory, reset positions and forces
                forces = _traj_forces[index];
                positions = _traj_positions[index];
                index++;

                //Record CV starting values
                for(size_t i = 0; i < _cv_start.size(); i++)
                {
                    _cv_start[i] = cvs[i]; 
                }
            }

            if((_iterator - _initialize_steps - _restrained_steps) % _swarm_length == 0)
            {
                //End of trajectory, harvest drift
                for(size_t i = 0; i < _cv_drift.size(); i++)
                {
                    _cv_drift[i] += cvs[i] - _cv_start[i]; //Add up drifts, average later
                }
            }
            _iterator++;
        }
        else
        {
            //Average drift
            for(size_t i = 0; i < _cv_drift.size(); i++)
            {
                _cv_drift[i] /= _number_trajectories;
            }

            //Evolve CVs, reparametrize, and reset vectors
            _currentiter++;

            StringUpdate();

            PrintString(cvs);

            _iterator == 0;

            for(size_t i = 0; i < _cv_drift.size(); i++)
            {
                _cv_drift[i] = 0; 
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
        //Write node, iteration, centers of the string and current CV value to output file
        _stringout.precision(8);
        _stringout << _mpiid << " " << _currentiter << " ";

        for(size_t i = 0; i < _centers.size(); i++)
        {
            _stringout << _centers[i] << " " << CV[i]->GetValue() << " ";
        }
        _stringout << std::endl;

        //Write same info to terminal, omit current CV value
        std::cout << _mpiid < " " << _currentiter << " ";
        for(size_t i = 0; i < _centers.size(); i++)
        {
            std::cout << _centers[i] << " ";
        }
        std::cout << std::endl;
    }

    void Swarm::StringUpdate()
    {
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
        std::vector<double> lower_cv_neighbor, upper_cv_neighbor;
        lower_cv_neighbor.resize(centerssize);
        upper_cv_neighbor.resize(centerssize);

        //Set up nodes to receive from their backward neighbor and send to their forward neighbor, wrapping around at the string ends
        if(_mpiid == 0)
        {
            sendneighbor = 1;
            recvneighbor = _world.size() - 1;
        }
        else if (_mpiid = _world.size() - 1)
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

        cvs_new.resize(centersize);
        alpha_star_vector.resize(_numnodes);
        cvs_new_vector.resize(_numnodes);

        //Evolve CVs with average drift
        for(i = 0; i < cvs_new.size(); i++)
        {
            cvs_new[i] = _centers[i] + _cv_drift[i];
        }

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

        for(i = 0; i < centerssize, i++)
        {
            //Take from new CV values in one dimension into a vector
            //cvs_new_vector is thus the CV value in one dimension at each node
            mpi::all_gather(_world, cvs_new[i], cvs_new_vector);
            //spl.setpoints(alpha_star_vector, vector of all new points in one particular dimension)
            spl.set_points(alpha_star_vector, cvs_newv); //Spline initialized with points on uneven mesh
            //Update CV values to lie on a regular mesh
            _centers[i] = spl(_alpha); 
        }
    }
}
