/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Jonathan K. Whitmer <jwhitme1@nd.edu>
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
#include "ElasticBand.h"
#include "CVs/CVManager.h"
#include "Snapshot.h" 
#include <math.h>
#include <iostream>

namespace SSAGES
{
	// Post-integration hook.
	void ElasticBand::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		auto& forces = snapshot->GetForces();
		auto cvs = cvmanager.GetCVs(cvmask_);

		// Apply umbrella to cvs
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			// Get current CV and gradient.
			auto& cv = cvs[i];
			auto& grad = cv->GetGradient();

			// Compute dV/dCV.
			auto diff = cv->GetDifference(centers_[i]);
			auto D = cvspring_[i] * diff;

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
				for(int k = 0; k < forces[j].size(); ++k)
					forces[j][k] -= D*grad[j][k];

			// If not equilibrating and has evolved enough steps,
			// generate the gradient
			if(iterator_ >= equilibrate_ && iterator_ % evolution_ == 0)
			{
			  //The force should include the cvspring; didn't before.
				newcenters_[i] += D;
				nsampled_++;
			}
		}

		// Restart iteration and zero gradients when moving onto
		// next elastic band iteration
		if(iterator_ > (equilibrate_ + evolution_ * nsamples_))
		{	
		       StringUpdate();
		       CheckEnd(cvs);
		       UpdateWorldString(cvs); 
		       PrintString(cvs); 

			iterator_ = 0;

			for(size_t ii = 0; ii < newcenters_.size(); ii++)
				newcenters_[ii] = 0;

			nsampled_ = 0;
			iteration_++;
		}

		iterator_++;
	}

	void ElasticBand::StringUpdate()
	{
                double dot=0, norm=0, dist1 = 0, dist2 = 0, springfrc = 0;
		std::vector<double> lcv0, ucv0, tngt;
		lcv0.resize(centers_.size(), 0);
		ucv0.resize(centers_.size(), 0);

		GatherNeighbors(&lcv0, &ucv0);

		tngt.resize(centers_.size());

		for(size_t ii = 0; ii<centers_.size(); ii++)
		{
			tngt[ii] = ucv0[ii] - lcv0[ii];
			norm+=tngt[ii]*tngt[ii];
		}

		norm=sqrt(norm);

		for(size_t ii = 0; ii < centers_.size(); ii++)  {
			tngt[ii] /= norm;
			newcenters_[ii] /= ((double) nsampled_);
			dot+=newcenters_[ii]*tngt[ii];

			//get two distances
			dist1 = (ucv0[ii]-centers_[ii])*(ucv0[ii]-centers_[ii]));
			dist2 = (centers_[ii] - ucv0[ii])*(centers_[ii] - ucv0[ii]);
      		}

  dist1 = sqrt(dist1);
  dist2 = sqrt(dist2);

		// Evolution of the images within the band
		// Endpoints evolve due to gradient alone.
		for(size_t ii = 0; ii < centers_.size(); ii++)
		{
			if((mpiid_ != 0) && (mpiid_ != world_.size()-1))
			{
			        //subtract out tangent from "real" force
				newcenters_[ii] -= dot*tngt[ii];
				
				//centers evolve according to perpendicular "real" force and parallel "spring" force
				centers_[ii]     = centers_[ii] + tau_ * (newcenters_[ii] + stringspring_ * (dist1 - dist2) * tngt[ii]);
			}
			else
			{
				centers_[ii] = centers_[ii] + tau_ * newcenters_[ii];
			}
		}
	}
}
