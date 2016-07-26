#include "Umbrella.h"

#include <iostream>

namespace SSAGES
{
	//! Value of the harmonic potential. Helper function.
	/*!
	 * \param k Spring constant.
	 * \param x0 Equilibrium extension.
	 * \param x Current extension.
	 * \return Harmonic potential at extension x.
	 */
	double spring(double k, double x0, double x)
	{
		return 0.5 * k * (x - x0) * (x - x0);
	}

	void Umbrella::PreSimulation(Snapshot*, const CVList& cvs)
	{
		if(_comm.rank() ==0)
		{
			char file[1024];
			sprintf(file, "node-%d.log", _world.rank());
		 	_umbrella.open(file);
		 	_currentiter = 0;
		 }
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
			auto D = _kspring[i]*(cv->GetDifference(_centers[i]));

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
				for(size_t k = 0; k < forces[j].size(); ++k)
					forces[j][k] -= D*grad[j][k];
		}
		PrintUmbrella(cvs);
		_currentiter++;
	}

	void Umbrella::PostSimulation(Snapshot*, const CVList&)
	{
		if(_comm.rank() ==0)
		{
			_umbrella.close();
		}
	}

	void Umbrella::PrintUmbrella(const CVList& CV)
	{
		if(_comm.rank() ==0)
		{
			_umbrella.precision(8);
			_umbrella << _currentiter << " ";

			for(size_t jj = 0; jj < _centers.size(); jj++)
				_umbrella<< _centers[jj] << " " << CV[jj]->GetValue()<< " "; 

			_umbrella<<std::endl;
		}
	}
}