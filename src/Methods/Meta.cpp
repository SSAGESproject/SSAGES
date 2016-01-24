#include "Meta.h"
#include <math.h>
#include <iostream>

namespace SSAGES
{
	double gaussian(double dx, double sigma)
	{
		double arg = (dx * dx) / (2. * sigma * sigma);
		return exp(-arg);
	}

	double gaussianDerv(double dx, double sigma)
	{
		double arg =  (dx * dx) / (2. * sigma * sigma);
		double pre = - dx / (sigma * sigma);
		return pre * exp(-arg);
	}

	void Meta::PreSimulation(Snapshot*, const CVList& cvs)
	{
		// TODO: Check widths against CV list size.

		// Get initial values of collective variables.
		for(auto& cv : cvs)
			_cvs.push_back(cv->GetValue());
	}

	void Meta::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		// Add hills when needed.
		if(snapshot->GetIteration() % _hillfreq == 0)
			AddHill(cvs);

		// Always calculate the current bias.
		CalcBiasForce();

		// TODO: Bounds check needs to go in somewhere.

		// Take each CV and add its biased forces to the atoms
		// using the chain rule.
		auto& forces = snapshot->GetForces();
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			auto& grad = cvs[i]->GetGradient();

			// Update the forces in snapshot by adding in the force bias from each
			// CV to each atom based on the gradient of the CV.
			for (size_t j = 0; j < forces.size(); ++j)
				for(size_t k = 0; k < forces[j].size(); ++k)
					forces[j][k] -= _derivatives[i]*grad[j][k];
		}
	}

	void Meta::AddHill(const CVList& cvs)
	{
		// Get cv values.
		for(size_t i = 0; i < cvs.size(); ++i)
			_cvs[i] = cvs[i]->GetValue();

		// Note: emplace_back constructs a hill in-place.
		_hills.emplace_back(_cvs, _widths, _height);
	}

	void Meta::CalcBiasForce()
	{	
		// Reset bias and derivatives.
		_bias = 0;
		for(size_t i = 0; i < _derivatives.size(); ++i)
			_derivatives[i] = 0;

		// Initialize vectors and reserve memory for calculation.
		std::vector<double> tder, dx; 
		tder.reserve(_cvs.size());
		dx.reserve(_cvs.size());

		// Loop through hills and calculate the bias force.
		for(auto& hill : _hills)
		{
			auto n = hill.center.size();
			auto tbias = 1.;
			
			// Resize vectors. 
			tder.resize(n, 1.0);
			dx.resize(n, 0);
			
			// Initialize dx and tbias.
			for(size_t i = 0; i < n; ++i)
			{
				dx[i] = _cvs[i] - hill.center[i];
				tbias *= gaussian(dx[i], hill.width[i]);
			}

			for(size_t i = 0; i < n; ++i)
				for(size_t j = 0; j < n; ++j)
				{
					if(j != i) 
						tder[i] *= gaussian(dx[j], hill.width[j]);
					else
						tder[i] *= gaussianDerv(dx[j], hill.width[j]);
				}

			_bias += _height * tbias;
			for(size_t i = 0; i < n; ++i)
				_derivatives[i] += _height*tder[i];
		}
	}
}