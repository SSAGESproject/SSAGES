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
#include <math.h>
#include <iostream>
#include <algorithm>
#include "../spline.h"

namespace SSAGES
{ 

	// Check whether CV values are within their respective Voronoi cell in CV space
	bool FiniteTempString::InCell(const CVList& cvs) const
	{
		std::vector<double> dists (_numnodes, 0);

		// Record the difference between all cvs and all nodes
		for (int i = 0; i < _numnodes; i++)
			for(size_t j = 0; j < cvs.size(); j++)
				dists[i]+= pow(cvs[j]->GetDifference(_worldstring[i][j]),2);
		
		if(std::min_element(dists.begin(), dists.end()) - dists.begin() == _mpiid)
			return true;

		return false;
	}

	// Post-integration hook.
	void FiniteTempString::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
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
        	auto& IDs = snapshot->GetAtomIDs();

        	for(size_t i = 0; i < _prev_IDs[0].size(); i++)
        	{
        		auto localindex = snapshot->GetLocalIndex(_prev_IDs[0][i]);
        		if(localindex!= -1)
        		{
        			Pos[localindex][0] = _prev_positions[0][i*3];
        			Pos[localindex][1] = _prev_positions[0][i*3 + 1];
        			Pos[localindex][2] = _prev_positions[0][i*3 + 2];
        		}
        	}                       
        }
        for(auto& cv : cvs)
        {
            //Trigger a rebuild of the CVs since we reset the positions
            cv->Evaluate(*snapshot);
        }
        insidecell = InCell(cvs); 

        MPI_Allreduce(MPI::IN_PLACE, &_run_umbrella, 1, MPI::BOOL, MPI::LOR, _world); 
		if(_run_umbrella)
		{ 
			if(_umbrella_iter == _min_num_umbrella_steps)
			{
				if(insidecell)
                {
                    _run_umbrella = false;
                    reset_for_umbrella = true;

                    //This node is done initializing; so store this snapshot
                    _prev_positions[0] = snapshot->SerializePositions();
                    _prev_IDs[0] = snapshot->SerializeIDs();
                }
				_umbrella_iter = 1;	
			}
			else
			{
                if(!reset_for_umbrella)
                {
                    for(size_t i = 0; i < cvs.size(); i++)
                    {
                        // Get current cv and gradient
                        auto& cv = cvs[i];
                        auto& grad = cv->GetGradient();

                        // Compute dV/dCV
                        auto D = _cvspring[i]*(cv->GetDifference(_centers[i]));

                        // Update forces
                        for(size_t j = 0; j < forces.size(); j++)
                                forces[j] -= D*grad[j];
                    }
                    _umbrella_iter++;    
                }
                else
                {
                    _run_umbrella = false;
                }
			}
            return;
		}
        else
        {
            reset_for_umbrella = false;
        }

		if(!insidecell)
		{
			for(auto& force : forces)
				force.setZero();

			auto& Pos = snapshot->GetPositions();
			auto& IDs = snapshot->GetAtomIDs();

			for(size_t i = 0; i < _prev_IDs[0].size(); i++)
			{
				auto localindex = snapshot->GetLocalIndex(_prev_IDs[0][i]);
				if(localindex!= -1)
				{
					Pos[localindex][0] = _prev_positions[0][i*3];
					Pos[localindex][1] = _prev_positions[0][i*3 + 1];
					Pos[localindex][2] = _prev_positions[0][i*3 + 2];
				}
			} 
			
			// Calculate running averages for each CV at each node based on previous CV
			for(size_t i = 0; i < _newcenters.size(); i++)
			{
				_newcenters[i] = _newcenters[i] * (_iteration * _blockiterations + _iterator - 1) + _prev_CVs[i];
				_newcenters[i] /= (_iteration * _blockiterations + _iterator);
			}
		}
		else
		{
			// Calculate running averages for each CV at each node 
			for(size_t i = 0; i < _newcenters.size(); i++)
			{
				_newcenters[i] = _newcenters[i] * (_iteration * _blockiterations + _iterator - 1) + cvs[i]->GetMinimumImage(_centers[i]);
				_newcenters[i] /= (_iteration * _blockiterations + _iterator);
			}

            _prev_CVs.clear();
            for(size_t i = 0; i < _centers.size(); i++)
            {
                _prev_CVs.push_back(cvs[i]->GetMinimumImage(_centers[i]));
            }
			_prev_positions[0] = snapshot->SerializePositions();
			_prev_IDs[0] = snapshot->SerializeIDs();
		}

		// Update the string, every _blockiterations string method iterations
		if(_iterator % _blockiterations == 0)
		{
            MPI_Barrier(_world);
	        StringUpdate();
            CheckEnd(cvs);
            MPI_Barrier(_world);
			UpdateWorldString(cvs); 
            PrintString(cvs);

			_iterator = 0;
			_iteration++;

			if(!InCell(cvs))
            {
                _run_umbrella = true;
                reset_for_umbrella = false; 
            }
			MPI_Allreduce(MPI::IN_PLACE, &_run_umbrella, 1, MPI::BOOL, MPI::LOR, _world);
		}

		_iterator++;
	}

	void FiniteTempString::StringUpdate()
	{
		std::vector<double> lcv0, ucv0;
		lcv0.resize(_centers.size(), 0);
		ucv0.resize(_centers.size(), 0);

		GatherNeighbors(&lcv0, &ucv0);
		
		// Update node locations toward running averages:
		for(size_t i = 0; i < _centers.size(); i++)
		{
			if(_mpiid == 0 || _mpiid == _numnodes - 1)
				_centers[i] =_centers[i] - _tau * (_centers[i] - _newcenters[i]);
			else
				_centers[i] = _centers[i] - _tau * (_centers[i] - _newcenters[i]) + 
					(_kappa * _numnodes * _tau * (ucv0[i] + lcv0[i] - 2 * _centers[i]));
		}

		GatherNeighbors(&lcv0, &ucv0);
		double alphastar = distance(_centers, lcv0);
		StringReparam(alphastar);
	}
}
