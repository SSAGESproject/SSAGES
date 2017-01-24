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
        double threshold = 0.2;
        const double eps = 0.0000000001;
        double diff;

        //On the first iteration, check that the CVs are within (threshold*100)% of the center value they're associated to
        for(size_t i = 0; i < cvs.size(); i++)
        {
            if(_centers[i] <= eps)
            {//e.g. if _centers[i] = 0
                diff = std::abs((cvs[i]->GetValue() - (_centers[i]+eps)) / ((cvs[i]->GetValue() + (eps + _centers[i]))/2.0));
            }
            else
            {
                diff = std::abs((cvs[i]->GetValue() - (_centers[i])) / ((cvs[i]->GetValue() + (_centers[i]))/2.0));
            }
            if(diff >= threshold)
            {
                return false; //proceed to initialize again
            }
        } 
        return true; //OK to move on to regular sampling
    }

    void Swarm::PostIntegration(Snapshot* snapshot, const CVList& cvs)
    {
         auto& forces = snapshot->GetForces();
         auto& positions = snapshot->GetPositions();
         auto& velocities = snapshot->GetVelocities();
         auto& atomids = snapshot->GetAtomIDs();

        if(snapshot_stored)
        {
            initialized = true;
        }
        else
        {
            initialized = CVInitialized(cvs); //Whether to initialize or not
        }
        
        original_initialized = initialized;

        if(initialized && !sampling_started)
        {
            if(!snapshot_stored)
            {
                _prev_positions[_index] = snapshot->SerializePositions();
                _prev_velocities[_index] = snapshot->SerializeVelocities();
                _prev_IDs[_index] = snapshot->SerializeIDs();
                snapshot_stored = true;
            }
        }
        MPI_Allreduce(MPI::IN_PLACE, &initialized, 1, MPI::BOOL, MPI::LAND, _world);
        if(!initialized && !sampling_started)
        {//Ensure CVs are initialized well
            //Do restrained sampling, and do not harvest trajectories 
            if(!original_initialized)
            {
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
                        for(int k = 0; k < forces[j].size(); k++)
                        {
                            forces[j][k] -= (double)D*grad[j][k];
                        }
                    }
                }
            }
            else
            {
                //Reset positions and forces, keeping them at their initialized value
                _index = 0;
                for(auto& force: forces)
                    force.setZero();

                for(size_t i = 0; i < _prev_IDs[_index].size(); i++)
                {
                    auto localindex = snapshot->GetLocalIndex(_prev_IDs[_index][i]);
                    if(localindex!= -1)
                    {
                        positions[localindex][0] = _prev_positions[_index][i*3];
                        positions[localindex][1] = _prev_positions[_index][i*3 + 1];
                        positions[localindex][2] = _prev_positions[_index][i*3 + 2];

                        velocities[localindex][0] = _prev_positions[_index][i*3];
                        velocities[localindex][1] = _prev_positions[_index][i*3 + 1];
                        velocities[localindex][2] = _prev_positions[_index][i*3 + 2];

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
                        for(int k = 0; k < forces[j].size(); k++)
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
                        _prev_positions[_index] = snapshot->SerializePositions();
                        _prev_velocities[_index] = snapshot->SerializeVelocities();
                        _prev_IDs[_index] = snapshot->SerializeIDs();
                        _index++;
                    }
                }
                if(_iterator == _initialize_steps + _restrained_steps)
                {
                    //Reset positions and forces before first call to unrestrained sampling
                    _index = 0;
                    for(auto& force: forces)
                        force.setZero();

                    for(size_t i = 0; i < _prev_IDs[_index].size(); i++)
                    {
                        auto localindex = snapshot->GetLocalIndex(_prev_IDs[_index][i]);
                        if(localindex!= -1)
                        {
                            positions[localindex][0] = _prev_positions[_index][i*3];
                            positions[localindex][1] = _prev_positions[_index][i*3 + 1];
                            positions[localindex][2] = _prev_positions[_index][i*3 + 2];

                            velocities[localindex][0] = _prev_positions[_index][i*3];
                            velocities[localindex][1] = _prev_positions[_index][i*3 + 1];
                            velocities[localindex][2] = _prev_positions[_index][i*3 + 2];

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
                        _cv_drift[i] = (_cv_drift[i]*_index + cvs[i]->GetMinimumImage(_centers[i])  - _centers[i]) / (_index+1); //Calculate running average of drifts 
                    }
                    //Set up for next trajectory
                    _index++;

                    if(_index < _number_trajectories)
                    {
                        for(auto& force: forces)
                            force.setZero();

                        for(size_t i = 0; i < _prev_IDs[_index].size(); i++)
                        {
                            auto localindex = snapshot->GetLocalIndex(_prev_IDs[_index][i]);
                            if(localindex!= -1)
                            {
                                positions[localindex][0] = _prev_positions[_index][i*3];
                                positions[localindex][1] = _prev_positions[_index][i*3 + 1];
                                positions[localindex][2] = _prev_positions[_index][i*3 + 2];

                                velocities[localindex][0] = _prev_positions[_index][i*3];
                                velocities[localindex][1] = _prev_positions[_index][i*3 + 1];
                                velocities[localindex][2] = _prev_positions[_index][i*3 + 2];

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
                
                MPI_Barrier(_world);
                StringUpdate();
                CheckEnd(cvs);
                MPI_Barrier(_world);
			    UpdateWorldString(cvs); 
                PrintString(cvs);

                _iterator = 0;
                _index = 0;
                snapshot_stored = false;

                for(size_t i = 0; i < _cv_drift.size(); i++)
                {
                    _cv_drift[i] = 0; 
                }
            }
        } 
    }

    void Swarm::StringUpdate()
	{
		// Update node locations toward running averages:
		for(size_t i = 0; i < _centers.size(); i++)
		{
         
            _centers[i] = _centers[i] + _cv_drift[i];
        }
        
		std::vector<double> lcv0, ucv0;
		lcv0.resize(_centers.size(), 0);
		ucv0.resize(_centers.size(), 0);
     
		GatherNeighbors(&lcv0, &ucv0);
    
		double alphastar = distance(_centers, lcv0);
  
		StringReparam(alphastar);
  
	}
}

