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
		for (size_t i = 0; i < _numnodes; i++)
			for(size_t j = 0; j < cvs.size(); j++)
				dists[i]+=cvs[j]->GetDifference(_worldstring[i][j])*cvs[j]->GetDifference(_worldstring[i][j]);
		
		if(std::min_element(dists.begin(), dists.end()) - dists.begin() == _mpiid)
			return true;

		return false;
	}

	// Post-integration hook.
	void FiniteTempString::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		auto& forces = snapshot->GetForces();

		auto insidecell = InCell(cvs);

		if(_run_umbrella)
		{
			if(insidecell && _umbrella_iter >= _min_num_umbrella_steps)
			{
				_run_umbrella = false;
				_umbrella_iter = 0;
			}
			else
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
				return;
			}
		}

		if(!insidecell)
		{
			for(auto& force : forces)
				force.setZero();

			SetPos(snapshot);
			
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
				_newcenters[i] = _newcenters[i] * (_iteration * _blockiterations + _iterator - 1) + cvs[i]->GetValue();
				_newcenters[i] /= (_iteration * _blockiterations + _iterator);
			}

			_prev_CVs.clear();
			for(auto&cv : cvs)
				_prev_CVs.push_back(cv->GetValue());

			StoreSnapshot(snapshot);
		}

		// Update the string, every _blockiterations string method iterations
		if(_iterator % _blockiterations == 0)
		{
			PrintString(cvs);
	        StringUpdate();
	        CheckEnd();
			UpdateWorldString();

			_iterator = 0;
			_iteration++;

			if(!InCell(cvs))
				_run_umbrella = true;
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
		double alphastar = sqdist(_centers, lcv0);
		StringReparam(alphastar);
	}
}
