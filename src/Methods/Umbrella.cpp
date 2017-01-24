/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
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
#include "Umbrella.h"

#include <iostream>

namespace SSAGES
{
	void Umbrella::PreSimulation(Snapshot* /* snapshot */, const CVList& /* cvs */)
	{
		if(_comm.rank() == 0)
		 	_umbrella.open(_filename.c_str(), std::ofstream::out | std::ofstream::app);
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
			auto center = GetCurrentCenter(snapshot->GetIteration(), i);
			auto D = _kspring[i]*(cv->GetDifference(center));

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
				forces[j] -= D*grad[j];
		}

		_iteration++;
		if(_iteration % _logevery == 0)
			PrintUmbrella(cvs);
	}

	void Umbrella::PostSimulation(Snapshot*, const CVList&)
	{
		if(_comm.rank() ==0)
			_umbrella.close();
	}

	void Umbrella::PrintUmbrella(const CVList& cvs)
	{
		if(_comm.rank() ==0)
		{
			_umbrella.precision(8);
			_umbrella << _iteration << " ";

			for(size_t i = 0; i < cvs.size(); i++)
				_umbrella << GetCurrentCenter(_iteration, i) << " " << cvs[i]->GetValue() << " "; 

			_umbrella << std::endl;
		}
	}
}
