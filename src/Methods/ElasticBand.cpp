#include "ElasticBand.h"
#include <math.h>
#include <iostream>

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
	 		_curr_field[i] = 0;
	 	}

	 	_currentiter = 1;
	 	_restartiter = 1;
	 	_nsampled = 0;

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
			auto D = _kspring[i] * (cv->GetValue() - _centers[i]);

			// Update forces.
			for(size_t j = 0; j < forces.size(); ++j)
			{
				for(size_t k = 0; k < forces[j].size(); ++k)
				{
					forces[j][k] -= D*grad[j][k];

					// If not equilibrating and has evolved enough steps,
					// generate the gradient
					if(_restartiter >= _equilibrate && _restartiter % _evolution == 0)
					{
						_gradient[i] += D*grad[j][k];
						sampled = true;
					}
				}
			}
		}

		if(sampled)
			_nsampled++;

		// Restart iteration and zero gradients when moving onto
		// next elastic band iteration
		if(_restartiter > (_equilibrate + _evolution * _nsamples))
		{
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

			PrintString(cvs);
		}

		if(_currentiter == _iterations)
			exit(0);

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
			_stringout<< _centers[jj] << " " << CV[jj]->GetValue() << " ";

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
			dot+=_gradient[ii]*tngt[ii]/_nsampled;
			norm+=tngt[ii]*tngt[ii];
		}

		norm=sqrt(norm);
		dot/=norm;


		for(ii = 0; ii < _centers.size(); ii++) 
			tngt[ii] /= norm;

		// Evolution of the images and reparametrirized of the string
		for(ii = 0; ii < _centers.size(); ii++)
		{
			_gradient[ii] /= ((double) _nsampled);

			if((_mpiid != 0) && (_mpiid != _world.size()-1))
			{
				_gradient[ii] -= dot*tngt[ii];
				next_field[ii] = _curr_field[ii] + _timestep * 
				(_gradient[ii] + _stringspring * (ucv0[ii] + lcv0[ii] - 2 * _curr_field[ii])/_timestep);
			}
			else
			{
				next_field[ii] = _curr_field[ii] + _timestep * _gradient[ii];
			}

			_curr_field[ii] = next_field[ii];
		}
	}
}