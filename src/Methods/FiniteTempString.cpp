/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ashley Guo <azguo@uchicago.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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
#include "FiniteTempString.h"
#include "CVs/CVManager.h"
#include "Snapshot.h" 
#include "spline.h"
#include <math.h>
#include <iostream>
#include <algorithm>

namespace SSAGES
{ 
	// Check whether CV values are within their respective Voronoi cell in CV space
	bool FiniteTempString::InCell(const CVList& cvs) const
	{
		std::vector<double> dists (numnodes_, 0);

		// Record the difference between all cvs and all nodes
		for (int i = 0; i < numnodes_; i++)
			for(size_t j = 0; j < cvs.size(); j++)
				dists[i]+= pow(cvs[j]->GetDifference(worldstring_[i][j]),2);
		
		if(std::min_element(dists.begin(), dists.end()) - dists.begin() == mpiid_)
			return true;

		return false;
	}

	// Post-integration hook.
	void FiniteTempString::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		auto cvs = cvmanager.GetCVs(cvmask_);
		auto& forces = snapshot->GetForces();
        bool insidecell;

        if(reset_for_umbrella)
        {
            //If the system was not going to run the umbrella anymore, reset its position
            for(auto& force : forces)
            {
                force.setZero();
            }
        	
        	auto& Pos = snapshot->GetPositions();

        	for(size_t i = 0; i < prev_IDs_[0].size(); i++)
        	{
        		auto localindex = snapshot->GetLocalIndex(prev_IDs_[0][i]);
        		if(localindex!= -1)
        		{
        			Pos[localindex][0] = prev_positions_[0][i*3];
        			Pos[localindex][1] = prev_positions_[0][i*3 + 1];
        			Pos[localindex][2] = prev_positions_[0][i*3 + 2];
        		}
        	}                       
        }
        for(auto& cv : cvs)
        {
            //Trigger a rebuild of the CVs since we reset the positions
            cv->Evaluate(*snapshot);
        }
        insidecell = InCell(cvs); 

        MPI_Allreduce(MPI_IN_PLACE, &run_umbrella_, 1, MPI_INT, MPI_LOR, world_); //Compares all run_umbrellas on all processors - if any is true, all are set to true

		if(run_umbrella_)
		{ 
			if(umbrella_iter_ == min_num_umbrella_steps_)
			{
                //After reaching the minimum umbrella sampling checkpoint, evaluate whether to keep going (per individual image)
                if(insidecell)
                {
                    //Stop umbrella sampling if inside cell
                    run_umbrella_ = false;
                    reset_for_umbrella = true;

                    //This node is done initializing; so store this snapshot
                    prev_positions_[0] = snapshot->SerializePositions();
                    prev_IDs_[0] = snapshot->SerializeIDs();
                }
				umbrella_iter_ = 1;	
			}
			else
			{
                if(!reset_for_umbrella)
                {
                    //If not at the minimum step checkpoint AND not a system which was already done, add restraining force
                    for(size_t i = 0; i < cvs.size(); i++)
                    {
                        // Get current cv and gradient
                        auto& cv = cvs[i];

                        // Compute dV/dCV
                        auto D = cvspring_[i]*(cv->GetDifference(centers_[i]));

                        // Update forces
                        cv->ApplyBias(D,*snapshot);
                    }
                    umbrella_iter_++; //Progress toward checkpoint 
                }
                else
                {
                    //End umbrella sampling of this node if it was already initialized
                    run_umbrella_ = false;
                }
			}
            return; //Don't do any sampling after restraining
		}
        else
        {
            //If run umbrella was false (only possible if every node was done) then set reset to false for next string method iteration (only proceed past umbrella sampling if no images require further umbrella sampling)
            reset_for_umbrella = false;
        }

		if(!insidecell)
		{
            //If last step took the system outside the Voronoi cell, reset the position and zero the force to reset the trajectory 
			for(auto& force : forces)
				force.setZero();

			auto& Pos = snapshot->GetPositions();

			for(size_t i = 0; i < prev_IDs_[0].size(); i++)
			{
				auto localindex = snapshot->GetLocalIndex(prev_IDs_[0][i]);
				if(localindex!= -1)
				{
					Pos[localindex][0] = prev_positions_[0][i*3];
					Pos[localindex][1] = prev_positions_[0][i*3 + 1];
					Pos[localindex][2] = prev_positions_[0][i*3 + 2];
				}
			} 
			
			// Calculate running averages for each CV at each node based on previous CV
			for(size_t i = 0; i < newcenters_.size(); i++)
			{
				newcenters_[i] = newcenters_[i] * (iteration_ * blockiterations_ + iterator_ - 1) + prev_CVs_[i];
				newcenters_[i] /= (iteration_ * blockiterations_ + iterator_);
			}
		}
		else
		{
			// Calculate running averages for each CV at each node 
			for(size_t i = 0; i < newcenters_.size(); i++)
			{
				newcenters_[i] = newcenters_[i] * (iteration_ * blockiterations_ + iterator_ - 1) + cvs[i]->GetMinimumImage(centers_[i]);
				newcenters_[i] /= (iteration_ * blockiterations_ + iterator_);
			}
            
            //Store info about this step in prev_vectors
            prev_CVs_.clear();
            for(size_t i = 0; i < centers_.size(); i++)
            {
                prev_CVs_.push_back(cvs[i]->GetMinimumImage(centers_[i]));
            }
			prev_positions_[0] = snapshot->SerializePositions();
			prev_IDs_[0] = snapshot->SerializeIDs();
		}

		// Update the string, every blockiterations_ string method iterations
		if(iterator_ % blockiterations_ == 0)
		{
            MPI_Barrier(world_);
	        StringUpdate();
            CheckEnd(cvs);
            MPI_Barrier(world_);
			UpdateWorldString(cvs); 
            PrintString(cvs);

			iterator_ = 0; //Number of steps within current method iteration
			iteration_++; //Increment string method iteration

			if(!InCell(cvs))
            {
                run_umbrella_ = true;
                reset_for_umbrella = false; 
            }
            MPI_Allreduce(MPI_IN_PLACE, &run_umbrella_, 1, MPI_INT, MPI_LOR, world_); //Compares all run_umbrellas on all processors - if any is true, all are set to true	
		}

        iterator_++; //Iterate number of steps within current method iteration	
	}

	void FiniteTempString::StringUpdate()
	{
		std::vector<double> lcv0, ucv0;
		lcv0.resize(centers_.size(), 0);
		ucv0.resize(centers_.size(), 0);

		GatherNeighbors(&lcv0, &ucv0);
		
		// Update node locations toward running averages:
		for(size_t i = 0; i < centers_.size(); i++)
		{
			if(mpiid_ == 0 || mpiid_ == numnodes_ - 1)
				centers_[i] =centers_[i] - tau_ * (centers_[i] - newcenters_[i]);
			else
				centers_[i] = centers_[i] - tau_ * (centers_[i] - newcenters_[i]) + 
					(kappa_ * numnodes_ * tau_ * (ucv0[i] + lcv0[i] - 2 * centers_[i]));
		}

		GatherNeighbors(&lcv0, &ucv0);
		double alphastar = distance(centers_, lcv0);
		StringReparam(alphastar);
	}
}
