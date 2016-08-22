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
				dists[i]+=(cvs[j]->GetValue() - _worldstring[i][j])*(cvs[j]->GetValue() - _worldstring[i][j]);
		
		if(std::min_element(dists.begin(), dists.end()) - dists.begin() == _mpiid)
			return true;

		return false;
	}

	// Post-integration hook.
	void FiniteTempString::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		auto& forces = snapshot->GetForces();
		auto inside = InCell(cvs);

		if(!inside)
		{
			for(size_t i = 0; i < cvs.size(); i++)
			{
				// Get current cv and gradient
				auto& cv = cvs[i];
				auto& grad = cv->GetGradient();

				// Compute dV/dCV, slowly increase spring strength to ensure don't get stuck.
				auto D = _cvspring[i]*(1.0 + _spring_iter/100.0)*(cv->GetDifference(_centers[i]));

				// Update forces
				for(size_t j = 0; j < forces.size(); j++)
					for(size_t k = 0; k < forces[j].size(); k++)
						forces[j][k] -= D*grad[j][k];
			}
			_spring_iter++;
		}
		else
		{
			// Calculate running averages for each CV at each node 
			for(size_t i = 0; i < _newcenters.size(); i++)
			{
				_newcenters[i] = _newcenters[i] * (_iteration * _blockiterations + _iterator - 1) + cvs[i]->GetValue();
				_newcenters[i] /= (_iteration * _blockiterations + _iterator);
			}

			// Update the string, every _blockiterations string method iterations
			if(_iterator % _blockiterations == 0)
			{
				PrintString(cvs);
		        StringUpdate();
		        CheckEnd(cvs);
				UpdateWorldString();

				_iterator = 1;
				_iteration++;
			}
			else
			{
				_iterator++;
			}
			
			_spring_iter = 1;
		}
	}

	void FiniteTempString::StringUpdate()
	{
		std::vector<double> lcv0, ucv0;
		lcv0.resize(_centers.size(), 0);
		ucv0.resize(_centers.size(), 0);

		GatherNeighbors(&lcv0, &ucv0);

		double alphastar = sqdist(_centers, lcv0);

		// Update node locations toward running averages:
		for(size_t i = 0; i < _centers.size(); i++)
		{
			if(_mpiid == 0 || _mpiid == _numnodes - 1)
				_centers[i] =_centers[i] - _tau * (_centers[i] - _newcenters[i]);
			else
				_centers[i] = _centers[i] - _tau * (_centers[i] - _newcenters[i]) + 
					(_kappa * _numnodes * _tau * (ucv0[i] + lcv0[i] - 2 * _centers[i]));
		}

		StringReparam(alphastar);
	}
}
