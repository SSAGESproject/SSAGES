#pragma once

#include "Meta.h"
#include <math.h>

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

	void Meta::PreSimulation(Snapshot* snapshot, const CVList& cvs)
	{
		// Get initial values of collective variables.
		for(auto& cv : cvs)
			_cvs.push_back(cv->GetValue());
	}

	void Meta::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		// TODO : This is condition for adding new hills. Should be changed
		// to whatever is appropriate. If a user specified interval is needed
		// it can be passed into the constructor.
		if(snapshot->GetIteration() % 5 == 0)
			AddHill(cvs);

		// Always calculate the current bias force and chain rule.
		CalcBiasForce();
		ChainRule(cvs);

		// TODO: At this point we have computed forces (if I'm not mistaken).
		// To update them we go through our snapshot and update the forces 
		// on the atoms associated with the particular CVs. Below is just 
		// "rough" code. There is no private variable _forces. Make it if needed, 
		// and should also be filled out in ChainRule I guess. 
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			// Get references to current CV for readability.
			auto& cv = cvs[i];
			auto& ids = cv->GetAtomIDs();
			
			// Get reference to forces in snapshot for updating.
			auto& forces = snapshot->GetForces();

			// Go through atoms IDs for force update.
			for(size_t j = 0; j < ids.size(); ++j)
			{
				auto id = ids[j];
				// I'm pretend updating the forces on atom "j" for the current
				// CV, "i". "forces" is a vector of 3x1 arrays, where the id 
				// obtains from the ID list of the CV is the index of the manipulated
				// atom. 
				// In this case, I am pretending that _forces is a private variables of 
				// dimensions "cv count" x "cv atoms count" x 3. Note that I use index "j"
				// instead of "id" because the only entries in the _forces vector are those
				// corresponding to manipulated atoms.
				forces[id][0] += _forces[i][j][0];
				forces[id][1] += _forces[i][j][1];
				forces[id][2] += _forces[i][j][2];

				// NOTE: I of course do not knw of the forces are additive for each CV 
				// or if a net force is computed and done only once. In that case we do 
				// not need to loop through CV's. We need to maintain a list of atoms that
				// are manipulated by Meta (the union of all atoms manipulated by CVs) and
				// then just have a single loop going through forces and updating the
				// corresponding IDs. I can take care of this part if the rest is filled out. 
			}
			
		}
	}

	void Meta::AddHill(const CVList& cvs)
	{
		// Get cv values.
		for(size_t i = 0; i < cvs.size(); ++i)
			_cvs[i] = cvs[i]->GetValue();

		// TODO: Widths and height are so far uninitialized. Are these 
		// user-specified options passed in through the constructor? If so 
		// then it should be added. Are there multiple widths but 1 height? 
		// It seemed that way from the code so I put it in like that. Check Meta.h 
		// and adjust accordingly. 
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

	// Fill this in. I just demoed the methods but didn't do anything really. 
	void Meta::ChainRule(const CVList& cvs)
	{
		// The chain rule needs access to the atoms 
		// involved in the CV. These are available here:
		for(int i = 0; i < cvs.size(); ++i)
		{
			// We just take a reference to the current CV for readability.
			auto& cv = cvs[i];

			// This vector contains the CV gradient which should be of length 
			// the number of atoms contributing to the CV. 
			auto& grad = cv->GetGradient();

			// Boundaries are also available. This should be of length 2. 
			auto& bound = cv->GetBoundaries();

			// List of atom IDs contributing to CV. This should be 
			// the same length as grad.
			auto& ids = cv->GetAtomIDs();

			// Do we need access to the forces on the atoms at this point? or later?
			// If needed pass snapshot* as an argument to ChainRule. If not see 
			// Postintegration method for next step. 
		}
	}
}