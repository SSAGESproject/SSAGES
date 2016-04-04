#include "FiniteTempString.h"
#include <math.h>
#include <iostream>

namespace SSAGES
{ 
	// need helper function for calculating distances 
	double sqdist(vector<double>& x, vector<double>& y)
	{
		double distance = 0;
		for (size_t i = 0; i < x.size(); i++){
			distance += (x[i] - y[i]) * (x[i] - y[i]);	
		}
		
		return distance;
	}

	// Pre-simulation hook.
	void FiniteTempString::PreSimulation(Snapshot* snapshot, const CVList& cvs)
	{
		// Open file for writing.
		auto& forces = snapshot->GetForces();
		_mpiid = snapshot->GetNode();
		char file[1024];
		sprintf(file, "node-%04d.log",_mpiid);
	 	_stringout.open(file);
	 	_gradient.resize(cvs.size());
	 	_curr_field.resize(cvs.size());
	 	_prev_positions.resize(snapshot->GetPositions().size());

	 	_iterator = 0;

	 	for(auto& position : _prev_positions)
	 		for(auto& dim : position)
	 			dim = 0;

		// initialize running averages
		for(size_t i = 0; i< _cvs.size(); ++i){
			_runavgs[i] = 0;
			_curr_field[i] = 0;
			_cv_prev[i] = _centers[i];
		}

		for(size_t i = 0; i< _centers.size(); i++){
			_alpha[i] = i / (_centers.size());
		}

		for(ii = 0; ii < _centers.size(); ii++)
		{
			MPI_Allgather(&_centers[ii], _centers[ii].size(), MPI_Float, 
				&_worldstring[ii], _worldstring[ii].size(), MPI_Float, _world);
		}

		MPI_Allgather(&_runavgs, _runavgs.size(), MPI_Float,
               &_worldaverages, _runavgs.size(), MPI_Float,
               _world);

	}

	// Post-integration hook.
	void FiniteTempString::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{

		auto& iter = snapshot->GetIteration();
		if(iter < _equilibrate)
			continue;

		std::vector<double> dists;
		dists.resize(_worldstring.size())
		auto& forces = snapshot->GetForces();
		auto& positions = snapshot->GetPositions();


		// Record the difference between all cvs and all nodes
		for (size_t i = 0; i < _worldstring.size(); i++)
		{
			dists[i] = 0;
			for(size_t j = 0; i < _cvs.size(); j++)
				dists[i]+=(_cvs[j] - _worldstring[i][j])*(_cvs[j] - _worldstring[i][j]);
		}

		// Check to see if cv is still in the correct voronio cell, if not reverse the move
		// Hard voronoi walls.
		// Reverse move by moving everything back, and setting force to zero
		// This should also 'Raondomize' the velocity, given the velocity will be that of the
		// unaccepted move. Does this follow detailed balance?

		if(std::min_element(dists) != _mpiid)
		{
			for(auto& force : forces)
				for(auto& xyz : force)
					xyz = 0.0;

			for(size_t i = 0; i < positions.size(); i++)
				for(size_t j = 0; j < positions[i].size(); j++)
					positions[i][j] = _prev_positions[i][j];
		}
		else
		{		
			for(size_t i =0; i<_cvs.size(); i++)
				_cv_prev[i] = _cvs[i]->GetValue();
		}

		// calculate running averages
		for (size_t i = 0; i < _runavgs.size(); i++)
		{
			// calculate running average for each node
			_runavgs[i] = _runavgs[i] * _iterator + _cv_prev[i];
			_runavgs[i] /= (_iterator + 1);
		}

		//write output if appropriate
		if (iter % _printfreq) {
			writeString();
		}

		if(_iterator > _blockiterations)
		{

			//update the string and reparameterize 
			StringUpdate();

			_iterator = 0;
			for (auto &cvavg : _runavgs)
				cvavg = 0;

			//wait for communication
			MPI_Barrier(comm);

			// update the umbrellas/string with new values
			for (size_t ii = 0; ii < _centers.size(); ii++)
				_centers[ii] = _curr_field[ii];
		}

		for(size_t i =0; i<_cvs.size(); i++)
			_cv_prev[i] = _cvs[i]->GetValue();

	}

	// Post-simulation hook.
	void FiniteTempString::PostSimulation(Snapshot*, const CVList&)
	{
		_stringout.close();
	}

	void FiniteTempString::PrintString(const CVList& CV)
	{
		_stringout.precision(8);
		_stringout << _mpiid << " "<< _currentiter << " "

		for(size_t jj = 0; jj < _centers.size(); jj++)
			_stringout<< _centers[jj], << " " << CV[jj]->GetValue() <<std::endl;

	}

	void FiniteTempString::StringUpdate()
	{
		size_t ii,jj;
		int centersize = _centers.size();
		std::vector<double> alpha_star;
		std::vector<double> cvs_new;

		std::vector<double> lcv0, ucv0;
		lcv0.resize(centersize);
		ucv0.resize(centersize);

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

		MPI_Sendrecv(&_centers[0], centersize, MPI_DOUBLE, sendneighbor, 1234,
			&lcv0[0], centersize, MPI_DOUBLE, recvneighbor, 1234, 
			_world, &status);

		MPI_Sendrecv(&_centers[0], centersize, MPI_DOUBLE, recvneighbor, 4321,
		       &ucv0[0], centersize, MPI_DOUBLE, sendneighbor, 4321, 
		       _world, &status);

		MPI_Allgather(&_runavgs, _runavgs.size(), MPI_Float, 
			&_worldaverages, _runavgs.size(), MPI_Float, _world);


		cvs_new.resize(_cvs.size());
		for(jj = 0; jj < cvs_new.size(); jj++)
		{
			if(_mpiid == 0 || _mpiid == _centers.size() - 1)
				cvs_new[jj] = _centers[jj] - dtau * (_centers[jj] - _runavgs[jj]);
			else
				cvs_new[jj] = _wolrdstring[ii][jj] - dtau * (_centers[jj] - _runavgs[jj]) + 
					(kappa * cvs_new.size() * dtau * 
					(ucv0[jj] + lcv0[jj] - 2 * _centers[jj]));
		}

		alpha_star.push_back(0);
		for(ii = 1; ii < _centers.size(); ii++)
			alpha_star.push_back() = alpha_star[ii-1] + sqrt(sqdist(_wolrdstring[ii], _wolrdstring[ii-1]));

		for(ii = 0; ii < _centers.size(); ii++)
			alpha_star[ii] /= alpha_star[_centers.size() - 1];

		for(ii = 0; ii < cvs_new.size(); ii++)
		{
			tk::spline spl.set_points(alpha_star[ii], cvs_new[ii]);
			for(jj = 0; jj < cvs_new[ii].size(); jj++)
				_wolrdstring[jj][ii] = spl(_alpha[jj]);
		}
	}