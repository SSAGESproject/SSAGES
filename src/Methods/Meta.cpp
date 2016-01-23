#pragma once

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
		// TODO : This is condition for adding new hills. Should be changed
		// to whatever is appropriate. If a user specified interval is needed
		// it can be passed into the constructor.
		if(snapshot->GetIteration() % _hillfreq == 0)
			AddHill(cvs);

		// Always calculate the current bias.
		CalcBiasForce();

		// Take each CV and add its biased forces to the atoms
		// using the chain rule.
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			// We just take a reference to the current CV and derivative for readability.
			auto& cv = cvs[i];
			auto& derivative = _derivatives[i];

			// This vector contains the CV gradient which should be of length 
			// the number of atoms. Atoms that do not contribute have a gradient of zero. 
			auto& grad = cv->GetGradient();

			// Boundaries are also available. This should be of length 2.
			auto& bound = cv->GetBoundaries();

			auto& forces = snapshot->GetForces();

			//Sanity check:
			if(forces.size() != grad.size())
			{
				std::cerr 
				<< "Error - cannot calculate dot product of mismatched matrix!" << std::endl
				<< "# atoms snap shot : " << forces.size() << ". # atoms cv->GetGradient : " 
				<< grad.size() << std::endl;
				exit(-1);
			}

			// Update the forces in snapshot by adding in the force bias from each
			// CV to each atom based on the gradient of the CV.
			for (size_t atomn = 0; atomn < forces.size(); ++atomn)
				for(size_t fxyz = 0; fxyz < forces[atomn].size(); ++fxyz)
					forces[atomn][fxyz] += derivative*grad[atomn][fxyz];
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
				// TODO: is hill.width[i] what is meant by sigmas?
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