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
            if(centers_[i] <= eps)
            {//e.g. if centers_[i] = 0
                diff = std::abs((cvs[i]->GetValue() - (centers_[i]+eps)) / ((cvs[i]->GetValue() + (eps + centers_[i]))/2.0));
            }
            else
            {
                diff = std::abs((cvs[i]->GetValue() - (centers_[i])) / ((cvs[i]->GetValue() + (centers_[i]))/2.0));
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
                prev_positions_[index_] = snapshot->SerializePositions();
                prev_velocities_[index_] = snapshot->SerializeVelocities();
                prev_IDs_[index_] = snapshot->SerializeIDs();
                snapshot_stored = true;
            }
        }
        MPI_Allreduce(MPI::IN_PLACE, &initialized, 1, MPI::BOOL, MPI::LAND, world_);
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
                    auto D = cvspring_[i]*(cv->GetDifference(centers_[i]));

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
                index_ = 0;
                for(auto& force: forces)
                    force.setZero();

                for(size_t i = 0; i < prev_IDs_[index_].size(); i++)
                {
                    auto localindex = snapshot->GetLocalIndex(prev_IDs_[index_][i]);
                    if(localindex!= -1)
                    {
                        positions[localindex][0] = prev_positions_[index_][i*3];
                        positions[localindex][1] = prev_positions_[index_][i*3 + 1];
                        positions[localindex][2] = prev_positions_[index_][i*3 + 2];

                        velocities[localindex][0] = prev_positions_[index_][i*3];
                        velocities[localindex][1] = prev_positions_[index_][i*3 + 1];
                        velocities[localindex][2] = prev_positions_[index_][i*3 + 2];

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
            if(iterator_ <= initialize_steps_ + restrained_steps_)
            {
                //Do restrained sampling, and do not harvest trajectories
                for(size_t i = 0; i < cvs.size(); i++)
                {
                    if(iterator_ == 0)
                    {
                        index_ = 0; //Reset index when starting
                    }
                    //Get current CV and gradient
                    auto& cv = cvs[i];
                    auto& grad = cv->GetGradient();

                    //Compute dV/dCV
                    auto D = cvspring_[i]*(cv->GetDifference(centers_[i]));

                    //Update forces
                    for(size_t j = 0; j < forces.size(); j++)
                    {
                        for(int k = 0; k < forces[j].size(); k++)
                        {
                            forces[j][k] -= (double)D*grad[j][k]; 
                        }
                    }
                }
                if(iterator_ > initialize_steps_)
                {
                    //Harvest a trajectory every harvest_length_ steps
                    if(iterator_ % harvest_length_ == 0)
                    {
                        prev_positions_[index_] = snapshot->SerializePositions();
                        prev_velocities_[index_] = snapshot->SerializeVelocities();
                        prev_IDs_[index_] = snapshot->SerializeIDs();
                        index_++;
                    }
                }
                if(iterator_ == initialize_steps_ + restrained_steps_)
                {
                    //Reset positions and forces before first call to unrestrained sampling
                    index_ = 0;
                    for(auto& force: forces)
                        force.setZero();

                    for(size_t i = 0; i < prev_IDs_[index_].size(); i++)
                    {
                        auto localindex = snapshot->GetLocalIndex(prev_IDs_[index_][i]);
                        if(localindex!= -1)
                        {
                            positions[localindex][0] = prev_positions_[index_][i*3];
                            positions[localindex][1] = prev_positions_[index_][i*3 + 1];
                            positions[localindex][2] = prev_positions_[index_][i*3 + 2];

                            velocities[localindex][0] = prev_positions_[index_][i*3];
                            velocities[localindex][1] = prev_positions_[index_][i*3 + 1];
                            velocities[localindex][2] = prev_positions_[index_][i*3 + 2];

                        }
                    }
                }
                iterator_++;
            }
            else if(iterator_ <= initialize_steps_ + restrained_steps_ + unrestrained_steps_)
            {
                //Launch unrestrained trajectories
                if((iterator_ - initialize_steps_ - restrained_steps_) % swarm_length_ == 0)
                {
                    //End of trajectory, harvest drift
                    for(size_t i = 0; i < cv_drift_.size(); i++)
                    { 
                        cv_drift_[i] = (cv_drift_[i]*index_ + cvs[i]->GetMinimumImage(centers_[i])  - centers_[i]) / (index_+1); //Calculate running average of drifts 
                    }
                    //Set up for next trajectory
                    index_++;

                    if(index_ < number_trajectories_)
                    {
                        for(auto& force: forces)
                            force.setZero();

                        for(size_t i = 0; i < prev_IDs_[index_].size(); i++)
                        {
                            auto localindex = snapshot->GetLocalIndex(prev_IDs_[index_][i]);
                            if(localindex!= -1)
                            {
                                positions[localindex][0] = prev_positions_[index_][i*3];
                                positions[localindex][1] = prev_positions_[index_][i*3 + 1];
                                positions[localindex][2] = prev_positions_[index_][i*3 + 2];

                                velocities[localindex][0] = prev_positions_[index_][i*3];
                                velocities[localindex][1] = prev_positions_[index_][i*3 + 1];
                                velocities[localindex][2] = prev_positions_[index_][i*3 + 2];

                            }
                        }
                    }
                }
                iterator_++;
            }
            else
            {
                //Evolve CVs, reparametrize, and reset vectors
                iteration_++;
                
                MPI_Barrier(world_);
                StringUpdate();
                CheckEnd(cvs);
                MPI_Barrier(world_);
			    UpdateWorldString(cvs); 
                PrintString(cvs);

                iterator_ = 0;
                index_ = 0;
                snapshot_stored = false;

                for(size_t i = 0; i < cv_drift_.size(); i++)
                {
                    cv_drift_[i] = 0; 
                }
            }
        } 
    }

    void Swarm::StringUpdate()
	{
		// Update node locations toward running averages:
		for(size_t i = 0; i < centers_.size(); i++)
		{
         
            centers_[i] = centers_[i] + cv_drift_[i];
        }
        
		std::vector<double> lcv0, ucv0;
		lcv0.resize(centers_.size(), 0);
		ucv0.resize(centers_.size(), 0);
     
		GatherNeighbors(&lcv0, &ucv0);
    
		double alphastar = distance(centers_, lcv0);
  
		StringReparam(alphastar);
  
	}
}

