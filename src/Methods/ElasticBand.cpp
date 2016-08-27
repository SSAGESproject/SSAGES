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
				for(size_t k = 0; k < forces[j].size(); ++k)
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
			PrintString(cvs);
	        StringUpdate();
	        CheckEnd(cvs);
			UpdateWorldString(); 

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
