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
    //Helper function to check if CVs are initialized correctly
    bool Swarm::CVInitialized(const CVList& cvs)
    {
        double threshold = 0.05;
        const double epsilon = 0.000000001;
        double diff;

        //On the first iteration, check that the CVs are within (threshold*100)% of the center value they're associated to
        for(size_t i = 0; i < cvs.size(); i++)
        {
            if(_centers[i] <= epsilon)
            {//e.g. if _centers[i] = 0
                diff = std::abs((cvs[i]->GetValue() - (_centers[i]+0.01)) / ((cvs[i]->GetValue() + (0.01 + _centers[i]))/2.0));
            }
            else
            {
                diff = std::abs((cvs[i]->GetValue() - (_centers[i])) / ((cvs[i]->GetValue() + (_centers[i]))/2.0));
            }
            if(diff >= threshold)
            {
                return true; //e.g. proceed to initialize again
            }
        }
        return false; //e.g. OK to move on to regular sampling
    }

    void Swarm::PostIntegration(Snapshot* snapshot, const CVList& cvs)
    {
        auto& forces = snapshot->GetForces();
        auto& positions = snapshot->GetPositions();
        auto& velocities = snapshot->GetVelocities();

        bool initialize; //Whether to initialize or not

        if(!sampling_started) 
        {
            initialize = CVInitialized(cvs);
        }
        else
        {
            initialize = false;
        }
        if(initialize && !sampling_started)
        {//On first pass, make sure CVs are initialized well
            //Do restrained sampling, and do not harvest trajectories 
            for(size_t i = 0; i < cvs.size(); i++)
            {
                //Get current CV and gradient
                auto& cv = cvs[i];
                auto& grad = cv->GetGradient();

                //Compute dV/dCV
                auto D = _cvspring[i]*(cv->GetDifference(_centers[i]));

                //Update forces
                for(size_t j = 0; j < forces.size(); j++)
                {
                    for(size_t k = 0; k < forces[j].size(); k++)
                    {
                        forces[j][k] -= (double)D*grad[j][k];
                    }
                }
            }
        }
        else
        {
            if(!sampling_started)
            {
                sampling_started = true; //Flag to prevent unneeded umbrella sampling
            }
            if(_iterator <= _initialize_steps + _restrained_steps)
            {
                //Do restrained sampling, and do not harvest trajectories
                for(size_t i = 0; i < cvs.size(); i++)
                {
                    if(_iterator == 0)
                    {
                        _index = 0; //Reset index when starting
                    }
                    //Get current CV and gradient
                    auto& cv = cvs[i];
                    auto& grad = cv->GetGradient();

                    //Compute dV/dCV
                    auto D = _cvspring[i]*(cv->GetDifference(_centers[i]));

                    //Update forces
                    for(size_t j = 0; j < forces.size(); j++)
                    {
                        for(size_t k = 0; k < forces[j].size(); k++)
                        {
                            forces[j][k] -= (double)D*grad[j][k]; 
                        }
                    }
                }
                if(_iterator > _initialize_steps)
                {
                    //Harvest a trajectory every _harvest_length steps
                    if(_iterator % _harvest_length == 0)
                    {
                        for(size_t k = 0; k < positions.size(); k++)
                        {
                            for(size_t l = 0; l < positions[k].size(); l++)
                            {
                                _traj_positions[_index][k][l] = positions[k][l];
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
                        //Zero forces
                        for(auto& force: forces)
                            force.setZero();
                        //Then set positions
                        for(size_t k = 0; k < positions.size(); k++)
                        {
                            for(size_t l = 0; l < positions[k].size(); l++)
                            {
                                positions[k][l] = _traj_positions[_index][k][l];
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
                if((_iterator - _initialize_steps - _restrained_steps) % _swarm_length == 0)
                {
                    //End of trajectory, harvest drift
                    for(size_t i = 0; i < _cv_drift.size(); i++)
                    {
                        _cv_drift[i] = (_cv_drift[i]*_index + cvs[i]->GetValue()  - _centers[i]) / (_index+1); //Calculate running average of drifts
                    }
                    //Set up for next trajectory
                    _index++;
                    if(_index < _number_trajectories)
                    {
                        //Start of trajectory, reset positions and forces
                        {
                            //Zero forces
                            for(auto& force : forces)
                                force.setZero();
                            //Then set positions
                            for(size_t k = 0; k < positions.size(); k++)
                            {
                                for(size_t l = 0; l < positions[k].size(); l++)
                                {
                                    positions[k][l] = _traj_positions[_index][k][l];
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
            }
            else
            {
                //Evolve CVs, reparametrize, and reset vectors
                _iteration++;


                PrintString(cvs);
                StringUpdate();
                CheckEnd(cvs);
                UpdateWorldString();

                _iterator = 0;
                _index = 0;

                for(size_t i = 0; i < _cv_drift.size(); i++)
                {
                    _cv_drift[i] = 0; 
                }
            }
        }
    }

    void Swarm::StringUpdate()
	{
		std::vector<double> lcv0, ucv0;
		lcv0.resize(_centers.size(), 0);
		ucv0.resize(_centers.size(), 0);

		GatherNeighbors(&lcv0, &ucv0);

		double alphastar = sqdist(_centers, lcv0);

		// Update node locations toward running averages:
		for(size_t i = 0; i < _centers.size(); i++)
		{
            //std::cout << _centers[i] << " " << _cv_drift[i] << std::endl;
            _centers[i] = _centers[i] + _cv_drift[i];
        }
		StringReparam(alphastar);
	}
}

