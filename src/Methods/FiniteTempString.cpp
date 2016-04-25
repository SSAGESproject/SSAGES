#include "FiniteTempString.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include "../spline.h"

namespace mpi = boost::mpi;
namespace SSAGES
{ 
	// need helper function for calculating distances 
	double sqdist(std::vector<double>& x, std::vector<double>& y)
	{
		double distance = 0;
		for (size_t i = 0; i < x.size(); i++)
		{
			distance += (x[i] - y[i]) * (x[i] - y[i]);	
		}
		
		return distance;
	}

	bool FiniteTempString::InCell(const CVList& cvs)
	{
		std::vector<double> dists;
		dists.resize(_worldstring[0].size());
		// Record the difference between all cvs and all nodes
		for (size_t i = 0; i < _numnodes; i++)
		{
			dists[i] = 0;
			for(size_t j = 0; j < cvs.size(); j++)
				dists[i]+=(cvs[j]->GetValue() - _worldstring[j][i])*(cvs[j]->GetValue() - _worldstring[j][i]);
		}

		if(std::min_element(dists.begin(), dists.end()) - dists.begin() != _mpiid)
			return false;
		else
			return true;

	}

	// Pre-simulation hook.
	void FiniteTempString::PreSimulation(Snapshot* snapshot, const CVList& cvs)
	{
		// Open file for writing.
		_mpiid = snapshot->GetWalkerID();
		auto& positions = snapshot->GetPositions();
		char file[1024];
		sprintf(file, "node-%04d.log",_mpiid);
	 	_stringout.open(file);
	 	_prev_positions.resize(positions.size());
	 	_worldstring.resize(_centers.size());
	 	_runavgs.resize(_centers.size());
	 	_cv_prev.resize(_centers.size());
	 	_SMD_centers.resize(_centers.size());
	 	_SMD_lengths.resize(_centers.size());

	 	_iterator = 0;

	 	for(size_t i=0; i<positions.size();i++)
	 	{
	 		for(size_t j=0; j<positions[i].size();j++)
	 			_prev_positions[i][j] = positions[i][j];
	 	}

	 	_alpha = _mpiid / (_numnodes - 1.0);

		// initialize running averages
		for(size_t i = 0; i< _centers.size(); ++i){
			_worldstring[i].resize(_numnodes);
			_runavgs[i] = 0;
			_cv_prev[i] = cvs[i]->GetValue();

			// gathers into _worldstring, where it is cv followed by node
			mpi::all_gather(_world, _centers[i], _worldstring[i]);
		}

		_run_SMD = !InCell(cvs);
	}

	// Post-integration hook.
	void FiniteTempString::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{

		auto& forces = snapshot->GetForces();
		auto& positions = snapshot->GetPositions();


		// Check to see if cv is still in the correct voronio cell, if not reverse the move
		// Hard voronoi walls.
		// Reverse move by moving everything back, and setting force to zero
		// This should also 'Raondomize' the velocity, given the velocity will be that of the
		// unaccepted move. Does this follow detailed balance?
		if(_run_SMD)
		{
			double dist = 0;
			double length = 0;
			for(size_t i = 0; i < cvs.size(); ++i)
			{
				auto& cv = cvs[i];
				auto& center = _centers[i];

				dist += (center - cv->GetValue())*(center - cv->GetValue());
				length += _SMD_lengths[i]*_SMD_lengths[i]; 
			}

			if(sqrt(dist) < 2*sqrt(length))
				_run_SMD = false;
		}
		if(_run_SMD)
		{
			// Run some SMD code
			double spring = 10;

			// Set up vector for SMD
			if(_cv_inside_iterator == 0)
			{
				for(size_t i = 0; i < cvs.size(); ++i)
				{
					auto& cv = cvs[i];
					auto& center = _centers[i];

					_SMD_lengths[i] = (center - cv->GetValue())/2000.0;
					_SMD_centers[i] = cv->GetValue() + _SMD_lengths[i];
				}
			}
			else if(_cv_inside_iterator < 2000)
			{
				for(size_t i = 0; i < _SMD_centers.size(); ++i)
					_SMD_centers[i]+=_SMD_lengths[i];
			}
			else
			{
				spring += _cv_inside_iterator - 2000;
			}

			// Typical harmonic spring to adjust forces
			for(size_t i = 0; i < cvs.size(); ++i)
			{
				// Get current CV and gradient.
				auto& cv = cvs[i];
				auto& grad = cv->GetGradient();

				// Compute dV/dCV.
				auto D = spring*(cv->GetValue() - _SMD_centers[i]);

				// Update forces.
				for(size_t j = 0; j < forces.size(); ++j)
					for(size_t k = 0; k < forces[j].size(); ++k)
						forces[j][k] -= D*grad[j][k];
			}
			_cv_inside_iterator++;
		}
		else
		{
			auto inside = InCell(cvs);
			// Do everything else
			if(!inside)
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
				for(size_t i = 0; i < positions.size(); i++)
					for(size_t j = 0; j < positions[i].size(); j++)
						_prev_positions[i][j] = positions[i][j];

				for(size_t i =0; i < cvs.size(); i++)
					_cv_prev[i] = cvs[i]->GetValue();
			}

			// calculate running averages
			for (size_t i = 0; i < _runavgs.size(); i++)
			{
				// calculate running average for each node
				_runavgs[i] = _runavgs[i] * _iterator + _cv_prev[i];
				_runavgs[i] /= (_iterator + 1);
			}

			if(_iterator > _blockiterations)
			{

				// Write out the string to file
				PrintString(cvs);

				// Update the string and reparameterize 
				StringUpdate();

				_iterator = 0;
				for (auto &cvavg : _runavgs)
						cvavg = 0;

				for(size_t ii = 0; ii < _centers.size(); ii++)
					mpi::all_gather(_world, _centers[ii], _worldstring[ii]);

				_currentiter++;

				if(!InCell(cvs))
				{
					_run_SMD = true;
					_cv_inside_iterator = 0;
				}
			} 
			else 
			{
				_iterator++;
				_run_SMD = false;
			}
		}
	}

	// Post-simulation hook.
	void FiniteTempString::PostSimulation(Snapshot*, const CVList&)
	{
		_stringout.close();
	}

	void FiniteTempString::PrintString(const CVList& CV)
	{
		_stringout.precision(8);
		_stringout << _mpiid << " "<< _currentiter << " ";

		for(size_t jj = 0; jj < _centers.size(); jj++)
			_stringout<< _centers[jj] << " " << CV[jj]->GetValue()<< " "; 

		_stringout<<std::endl;

		std::cout << _mpiid << " "<< _currentiter << " ";

		for(size_t jj = 0; jj < _centers.size(); jj++)
			std::cout<< _centers[jj] << " " << _runavgs[jj]<< " "; 

		std::cout<<std::endl;

	}

	void FiniteTempString::StringUpdate()
	{
		size_t jj;
		int centersize = _centers.size();
		double alpha_star;
		std::vector<double> alpha_starv;
		std::vector<double> cvs_new;
		std::vector<double> cvs_newv; 

		std::vector<double> lcv0, ucv0;
		lcv0.resize(centersize);
		ucv0.resize(centersize);

		int sendneighbor, recvneighbor;
		MPI_Status status;

		if(_mpiid == 0)
		{
			sendneighbor = 1;
			recvneighbor = _world.size()-1;
		} 
		else if (_mpiid == _world.size()-1)
		{
			sendneighbor = 0;
			recvneighbor = _world.size() - 2;
		} 
		else 
		{
			sendneighbor = _mpiid + 1;
			recvneighbor = _mpiid - 1;
		}

		MPI_Sendrecv(&_centers[0], centersize, MPI_DOUBLE, sendneighbor, 1234,
			&lcv0[0], centersize, MPI_DOUBLE, recvneighbor, 1234, 
			_world, &status);

		MPI_Sendrecv(&_centers[0], centersize, MPI_DOUBLE, recvneighbor, 4321,
		       &ucv0[0], centersize, MPI_DOUBLE, sendneighbor, 4321, 
		       _world, &status);

		cvs_new.resize(centersize);
		alpha_starv.resize(_numnodes);
		cvs_newv.resize(_numnodes);

		for(jj = 0; jj < cvs_new.size(); jj++)
		{
			if(_mpiid == 0 || _mpiid == _numnodes - 1)
				cvs_new[jj] = _centers[jj] - _tau * (_centers[jj] - _runavgs[jj]);
			else
				cvs_new[jj] = _centers[jj] - _tau * (_centers[jj] - _runavgs[jj]) + 
					(_kappa * _numnodes * _tau * 
					(ucv0[jj] + lcv0[jj] - 2 * _centers[jj]));
		}

		if(_mpiid == 0)
			alpha_star = 0;
		else
			alpha_star = sqrt(sqdist(_centers, lcv0));

		mpi::all_gather(_world, alpha_star, alpha_starv);

		for(jj = 1; jj < alpha_starv.size(); jj++)
		{
			alpha_starv[jj] = alpha_starv[jj-1] + alpha_starv[jj];
		}

		for(jj = 1; jj < alpha_starv.size(); jj++){
			alpha_starv[jj] /= alpha_starv[_numnodes - 1];
		}

		tk::spline spl;

		for(jj = 0; jj < centersize; jj++)
		{
			mpi::all_gather(_world, cvs_new[jj], cvs_newv);
			//spl.setpoints(alpha_starv, vector of all new points **in one particular dimension**);
			spl.set_points(alpha_starv, cvs_newv);
			_centers[jj] = spl(_alpha);
		}
	}
}
