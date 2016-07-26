#include "Meta.h"
#include <math.h>
#include <iostream>

namespace SSAGES
{
	//! Evlauate Gaussian. Helper function.
	/*!
	 * \param dx x-value.
	 * \param sigma Width of Gaussian
	 * \return Value at x of Gaussian with center at zero and width sigma.
	 */
	double gaussian(double dx, double sigma)
	{
		double arg = (dx * dx) / (2. * sigma * sigma);
		return exp(-arg);
	}

	//! Evaluate Gaussian derivative. Helper function.
	/*!
	 * \param dx Value of x.
	 * \param sigma Width of Gaussian.
	 * \return Derivative at x of Gaussian with center at zero and width sigma.
	 */
	double gaussianDerv(double dx, double sigma)
	{
		double arg =  (dx * dx) / (2. * sigma * sigma);
		double pre = - dx / (sigma * sigma);
		return pre * exp(-arg);
	}

	// Pre-simulation hook.
	void Meta::PreSimulation(Snapshot*, const CVList& cvs)
	{
		// Open file for writing and allocate derivatives vector.
	 	_hillsout.open("hills.out");
		_derivatives.resize(cvs.size());	
	}

	// Post-integration hook.
	void Meta::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		// Add hills when needed.
		if(snapshot->GetIteration() % _hillfreq == 0)
			AddHill(cvs);

		// Always calculate the current bias.
		CalcBiasForce(cvs);

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
				for(size_t k = 0; k < 3; ++k)
					forces[j][k] -= _derivatives[i]*grad[j][k];
		}
	}

	// Post-simulation hook.
	void Meta::PostSimulation(Snapshot*, const CVList&)
	{
		_hillsout.close();
	}

	// Drop a new hill.
	void Meta::AddHill(const CVList& cvs)
	{
		std::vector<double> cvals;

		// Get CV values.
		for(size_t i = 0; i < cvs.size(); ++i)
			cvals.push_back(cvs[i]->GetValue());

		// Note: emplace_back constructs a hill in-place.
		_hills.emplace_back(cvals, _widths, _height);

		// Write hill to file.
		PrintHill(_hills.back());
	}

	//Ruthless pragmatism
	void Meta::PrintHill(const Hill& hill)
	{
		_hillsout.precision(8);
		for(auto& cv : hill.center)
			_hillsout << cv << " ";

		for(auto& w : hill.width)
			_hillsout << w << " ";

		_hillsout << _height << std::endl;
	}

	void Meta::CalcBiasForce(const CVList& cvs)
	{	
		// Reset bias and derivatives.
		_bias = 0;
		for(size_t i = 0; i < _derivatives.size(); ++i)
			_derivatives[i] = 0;

		// Loop through hills and calculate the bias force.
		for(auto& hill : _hills)
		{
			auto n = hill.center.size();
			std::vector<double> tder(n, 1.0), dx(n, 1); 
			auto tbias = 1.;
			
			// Initialize dx and tbias.
			for(size_t i = 0; i < n; ++i)
			{
				dx[i] = cvs[i]->GetDifference(hill.center[i]);
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
