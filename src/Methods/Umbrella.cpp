#include "Umbrella.h"

#include <iostream>

namespace SSAGES
{
	// Value of the harmonic potential. Helper function.
	double spring(double k, double x0, double x)
	{
		return 0.5 * k * (x - x0) * (x - x0);
	}

	// Derivative of harmonic potential. Helper function.
	double springDer(double k, double x0, double x)
	{
		return k * (x - x0);
	}

	void Umbrella::PreSimulation(Snapshot*, const CVList&)
	{
	}

	void Umbrella::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		// Compute the forces on the atoms from the CV's using the chain 
		// rule.
		auto& forces = snapshot->GetForces();
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			// Get current CV and gradient.
			auto& cv = cvs[i];
			auto& grad = cv->GetGradient();

			// Compute dV/dCV.
			auto D = springDer(_kspring[i], _centers[i], cv->GetValue());

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
				for(size_t k = 0; k < forces[j].size(); ++k)
					forces[j][k] -= D*grad[j][k];
		}
	}

	void Umbrella::PostSimulation(Snapshot*, const CVList&)
	{
	}
}