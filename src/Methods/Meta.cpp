/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Jonathan K. Whitmer <jwhitme1@nd.edu>
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
#include "Meta.h"
#include <math.h>
#include <iostream>
#include "Drivers/DriverException.h"

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
	void Meta::PreSimulation(Snapshot* snapshot, const CVList& cvs)
	{
		// Open file for writing and allocate derivatives vector.
		if(_world.rank() == 0)
			_hillsout.open("hills.out");

		auto n = snapshot->GetTargetIterations();
		n = n ? n : 1e5; // Pre-allocate at least something.
	
		_hills.reserve(n+1);
		_widths.reserve(n+1);
		_derivatives.resize(cvs.size());
		_tder.resize(cvs.size());
		_dx.resize(cvs.size());
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
		if(_world.rank() == 0)
			_hillsout.close();	
	}

	// Drop a new hill.
	void Meta::AddHill(const CVList& cvs)
	{
		int n = cvs.size();

		// Assume we have the same number of procs per walker.
		int nwalkers = _world.size()/_comm.size();

		// We need to exchange CV values across the walkers 
		// and to each proc on a walker.	
		std::vector<double> cvals(n*nwalkers, 0);

		if(_comm.rank() == 0)
		{
			for(auto i = 0, j = _world.rank()/_comm.size()*n; i < n; ++i,++j)
				cvals[j] = cvs[i]->GetValue();
		}

		// Reduce across all processors and add hills.
		MPI_Allreduce(MPI_IN_PLACE, cvals.data(), n*nwalkers, MPI_DOUBLE, MPI_SUM, _world);
		
		for(int i = 0; i < n*nwalkers; i += n)
		{
			std::vector<double> cval(cvals.begin() + i, cvals.begin() + i + n);
			_hills.emplace_back(cval, _widths, _height);
			
			// Write hill to file.
			if(_world.rank() == 0)
				PrintHill(_hills.back());
		}
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
		double bias = 0.;
		auto n = cvs.size();

		// Reset vectors.
		std::fill(_derivatives.begin(), _derivatives.end(), 0);

		// Loop through hills and calculate the bias force.
		for(auto& hill : _hills)
		{		
			auto tbias = 1.;
			std::fill(_tder.begin(), _tder.end(), 1.0);
			std::fill(_dx.begin(), _dx.end(), 1.0);
			
			for(size_t i = 0; i < n; ++i)
			{
				_dx[i] = cvs[i]->GetDifference(hill.center[i]);
				tbias *= gaussian(_dx[i], hill.width[i]);
			}

			for(size_t i = 0; i < n; ++i)
				for(size_t j = 0; j < n; ++j)
				{
					if(j != i) 
						_tder[i] *= gaussian(_dx[j], hill.width[j]);
					else
						_tder[i] *= gaussianDerv(_dx[j], hill.width[j]);
				}

			bias += _height * tbias;
			for(size_t i = 0; i < n; ++i)
				_derivatives[i] += _height*_tder[i];
		}
	}
}
