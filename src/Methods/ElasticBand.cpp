/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Jonathan K. Whitmer <jwhitme1@nd.edu>
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
#include "ElasticBand.h"
#include <math.h>
#include <iostream>
#include "../Drivers/DriverException.h"

namespace SSAGES
{
	//Helper functions here:


	// Each version of this method on each node will have different values in the CVList? QUESTION

	// Pre-simulation hook.
	void ElasticBand::PreSimulation(Snapshot* snapshot, const CVList& cvs)
	{
		// Open file for writing.
		_mpiid = snapshot->GetWalkerID();
		char file[1024];
		sprintf(file, "node-%04d.log",_mpiid);
	 	_stringout.open(file);
	 	_gradient.resize(cvs.size());
	 	_curr_field.resize(cvs.size());

	 	for(size_t i = 0; i < cvs.size(); ++i)
	 	{
	 		_gradient[i] = 0;
	 		_curr_field[i] = _centers[i];
	 	}

	 	_currentiter = 1;
	 	_restartiter = 1;
	 	_nsampled = 0;

	 	if(_centers.size() != _kspring.size() || _centers.size() != cvs.size())
	 		throw BuildException({"Number centers, kspring, and cvs are not the same!"});

	}

	// Post-integration hook.
	void ElasticBand::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{

		auto& forces = snapshot->GetForces();

		bool sampled = false;

		// Apply umbrella to cvs
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			// Get current CV and gradient.
			auto& cv = cvs[i];
			auto& grad = cv->GetGradient();

			// Compute dV/dCV.
			auto diff = cv->GetValue() - _centers[i];
			auto D = _kspring[i] * diff;

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
				for(size_t k = 0; k < forces[j].size(); ++k)
					forces[j][k] -= D*grad[j][k];

			// If not equilibrating and has evolved enough steps,
			// generate the gradient
			if(_restartiter >= _equilibrate && _restartiter % _evolution == 0)
			{
				_gradient[i] += diff;
				sampled = true;
			}
		}

		if(sampled)
			_nsampled++;

		// Restart iteration and zero gradients when moving onto
		// next elastic band iteration
		if(_restartiter > (_equilibrate + _evolution * _nsamples))
		{
			PrintString(cvs);

			//obtain the unparameterized string
			StringUpdate();

			//wait for communication
			MPI_Barrier(_world);

			// update the umbrellas/string with new values
			for (size_t ii = 0; ii < _centers.size(); ii++)
			{
				_centers[ii] = _curr_field[ii];
			} 

			_restartiter = 1;

			for(size_t ii = 0; ii < _gradient.size(); ii++)
				_gradient[ii] = 0;

			_currentiter++;
			_nsampled = 0;
		}

		_restartiter++;
		// TODO: Bounds check needs to go in somewhere.
	}

	// Post-simulation hook.
	void ElasticBand::PostSimulation(Snapshot*, const CVList&)
	{
		_stringout.close();
	}

	void ElasticBand::PrintString(const CVList& CV)
	{
 		_stringout.precision(8);
		_stringout << _mpiid << " "<< _currentiter << " ";

		for(size_t jj = 0; jj < _centers.size(); jj++)
		{
			_stringout<< _centers[jj] << " " << CV[jj]->GetValue() << " ";
		}

		_stringout << std::endl;

	}

	void ElasticBand::StringUpdate()
	{
		size_t ii;
		double dot=0, norm=0;
		int centersize = _centers.size();
		std::vector<double> lcv0, ucv0, tngt, next_field;
		lcv0.resize(centersize);
		ucv0.resize(centersize);
		tngt.resize(centersize);
		next_field.resize(centersize);

		int sendneighbor, recvneighbor;
		MPI_Status status;

		//get tangent vector
		//TODO: set nnodes somewhere
		if(_mpiid == 0){
			sendneighbor = 1;
			recvneighbor = _world.size()-1;
		} else if (_mpiid == _world.size()-1){
			sendneighbor = 0;
			recvneighbor = _world.size() - 2;
		} else {
			sendneighbor = _mpiid + 1;
			recvneighbor = _mpiid - 1;
		}

		MPI_Sendrecv(&_centers[0], _centers.size(), MPI_DOUBLE, sendneighbor, 1234,
			&lcv0[0], centersize, MPI_DOUBLE, recvneighbor, 1234, 
			_world, &status);


		MPI_Sendrecv(&_centers[0], _centers.size(), MPI_DOUBLE, recvneighbor, 4321,
		       &ucv0[0], centersize, MPI_DOUBLE, sendneighbor, 4321, 
		       _world, &status);

		for(ii = 0; ii<_centers.size(); ii++)
		{
			tngt[ii] = ucv0[ii] - lcv0[ii];
			norm+=tngt[ii]*tngt[ii];
		}

		norm=sqrt(norm);

		for(ii = 0; ii < _centers.size(); ii++)  {
			tngt[ii] /= norm;
		    _gradient[ii] /= ((double) _nsampled);
			dot+=_gradient[ii]*tngt[ii];
		}
		
		// Evolution of the images and reparametrirized of the string
		for(ii = 0; ii < _centers.size(); ii++)
		{
			if((_mpiid != 0) && (_mpiid != _world.size()-1))
			{
				_gradient[ii] -= dot*tngt[ii];
				next_field[ii] = _curr_field[ii] + _timestep * 
				(_gradient[ii] + _stringspring * (ucv0[ii] + lcv0[ii] - 2 * _curr_field[ii]));
			}
			else
			{
				next_field[ii] = _curr_field[ii] + _timestep * _gradient[ii];
			}

			_curr_field[ii] = next_field[ii];
		}


	}
}
