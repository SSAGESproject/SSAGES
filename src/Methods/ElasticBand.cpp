/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
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
#include <math.h>
#include <iostream>
#include "../Drivers/DriverException.h"

namespace SSAGES
{
	// Post-integration hook.
	void ElasticBand::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{

		auto& forces = snapshot->GetForces();

		// Apply umbrella to cvs
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			// Get current CV and gradient.
			auto& cv = cvs[i];
			auto& grad = cv->GetGradient();

			// Compute dV/dCV.
			auto diff = cv->GetDifference(_centers[i]);
			auto D = _cvspring[i] * diff;

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
				for(int k = 0; k < forces[j].size(); ++k)
					forces[j][k] -= D*grad[j][k];

			// If not equilibrating and has evolved enough steps,
			// generate the gradient
			if(_iterator >= _equilibrate && _iterator % _evolution == 0)
			{
				_newcenters[i] += diff;
				_nsampled++;
			}
		}

		// Restart iteration and zero gradients when moving onto
		// next elastic band iteration
		if(_iterator > (_equilibrate + _evolution * _nsamples))
		{	
	        StringUpdate();
            CheckEnd(cvs);
			UpdateWorldString(cvs); 
            PrintString(cvs); 

			_iterator = 0;

			for(size_t ii = 0; ii < _newcenters.size(); ii++)
				_newcenters[ii] = 0;

			_nsampled = 0;
			_iteration++;
		}

		_iterator++;
	}

	void ElasticBand::StringUpdate()
	{
		double dot=0, norm=0;
		std::vector<double> lcv0, ucv0, tngt;
		lcv0.resize(_centers.size(), 0);
		ucv0.resize(_centers.size(), 0);

		GatherNeighbors(&lcv0, &ucv0);

		tngt.resize(_centers.size());

		for(size_t ii = 0; ii<_centers.size(); ii++)
		{
			tngt[ii] = ucv0[ii] - lcv0[ii];
			norm+=tngt[ii]*tngt[ii];
		}

		norm=sqrt(norm);

		for(size_t ii = 0; ii < _centers.size(); ii++)  {
			tngt[ii] /= norm;
		    _newcenters[ii] /= ((double) _nsampled);
			dot+=_newcenters[ii]*tngt[ii];
		}
		
		// Evolution of the images and reparametrirized of the string
		for(size_t ii = 0; ii < _centers.size(); ii++)
		{
			if((_mpiid != 0) && (_mpiid != _world.size()-1))
			{
				_newcenters[ii] -= dot*tngt[ii];
				_centers[ii] = _centers[ii] + _tau * 
				(_newcenters[ii] + _stringspring * (ucv0[ii] + lcv0[ii] - 2 * _centers[ii]));
			}
			else
			{
				_centers[ii] = _centers[ii] + _tau * _newcenters[ii];
			}
		}
	}
}
