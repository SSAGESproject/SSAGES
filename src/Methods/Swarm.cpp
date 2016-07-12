#include "Swarm.h"
#include "../spline.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace mpi = boost::mpi;
namespace SSAGES
{
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

        index = 0;
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
                    _traj_positions[index];
                    _traj_forces[index]; 
                    index++;
                }
            }
            _iterator++;
        }
        else if(_iterator <= _initialize_steps + _restrained_steps + _unrestrained_steps*_number_trajectories)
        {
            //Launch unrestrained trajectories
            if(_iterator == _initialize_steps + _restrained_steps + 1)
            {
                index = 0; //Reset index when starting unrestrained trajectories
            }
            if((iterator - _initialize_steps - restrained_steps) % unrestrained_steps == 1)
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

            if((_iterator - _initialize_steps - _restrained_steps) % _unrestrained_steps == 0)
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
                _cv_drift = 0; 
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
        _stringout.precision(8);
        _stringout << _mpiid << " " << _currentiter << " ";

        for(size_t i = 0; i < _centers.size(); i++)
        {
            _stringout << _centers[i] << " " << CV[i]->GetValue() << " ";
        }
        _stringout << std::endl;

        std::cout << _mpiid < " " << _currentiter << " ";
        for(size_t i = 0; i < _centers.size(); i++)
        {
            std::cout << _centers[i] << " ";
        }
        std::cout << std::endl;
    }

    void Swarm::StringUpdate()
    {

    }
}
