/**        
- * This file is part of        
- * SSAGES - Suite for Advanced Generalized Ensemble Simulations        
- *     
- * Copyright 2016 Cody Bezik <bezik@uchicago.edu>      
- *     
- * SSAGES is free software: you can redistribute it and/or modify      
- * it under the terms of the GNU General Public License as published by        
- * the Free Software Foundation, either version 3 of the License, or       
- * (at your option) any later version.     
- *     
- * SSAGES is distributed in the hope that it will be useful,       
- * but WITHOUT ANY WARRANTY; without even the implied warranty of      
- * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
- * GNU General Public License for more details.        
- *     
- * You should have received a copy of the GNU General Public License       
- * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.     
- */
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
        auto& atomids = snapshot->GetAtomIDs();

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
                        for(size_t k = 0; k < atomids.size(); k++)
                        {
                            int current_id = atomids[k] - 1; //Atom ids are indexed from 1
                            for(size_t l = 0; l < positions[current_id].size(); l++)
                            {
                                _traj_positions[_index][current_id][l] = positions[k][l];
                            }
                            for(size_t l = 0; l < velocities[k].size(); l++)
                            {
                                _traj_velocities[_index][current_id][l] = velocities[k][l];
                            }
                            _traj_atomids[_index][current_id] = current_id + 1;
                        }
                        _index++;
                    }
                }
                if(_iterator == _initialize_steps + _restrained_steps)
                {
                    //Reset positions and forces before first call to unrestrained sampling
                    _index = 0;
                    ResetTrajectories(positions, velocities, forces, atomids);
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
                        ResetTrajectories(positions, velocities, forces, atomids);
                    }
                }
                _iterator++;
            }
            else
            {
                //Evolve CVs, reparametrize, and reset vectors
                _iteration++;
 
                StringUpdate();
                PrintString(cvs);
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

    void Swarm::ResetTrajectories(std::vector<Vector3>& positions, std::vector<Vector3>& velocities, std::vector<Vector3>& forces, Label& atomids)
    {
        //First, check that the system's atom IDs haven't changed before the swap - for now, abort if they have
        /*for(size_t k = 0; k < atomids.size(); k++)
        {
            if(atomids[k] != _traj_atomids[_index][k])
            {
                std::cout << "Your atom list has been made out of order!  The system cannot do a positional reset and thus must exit." << std::endl;
                _world.abort(-1);
            }
        }*/
        //Zero forces
        for(auto& force: forces)
            force.setZero();
        //Then set positions and velocities, taking care that the correct atom is actually reset
        for(size_t k = 0; k < atomids.size(); k++)
        {
            int current_id = atomids[k] - 1; //Atom ids are indexed from 1
            for(size_t l = 0; l < positions[k].size(); l++)
            {
                positions[k][l] = _traj_positions[_index][current_id][l];
            }
            for(size_t l = 0; l < velocities[k].size(); l++)
            {
                velocities[k][l] = _traj_velocities[_index][current_id][l];
            }
        }
    }
}

