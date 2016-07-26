#include "FiniteTempString.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include "../spline.h"

namespace mpi = boost::mpi;
namespace SSAGES
{ 
	//! Helper function for calculating distances
	/*!
	 * \param x List of coordinates.
	 * \param y List of coordinates.
	 * \return Sum of distances between the x-values and y-values.
	 */
	double sqdist(std::vector<double>& x, std::vector<double>& y)
	{
		double distance = 0;
		for (size_t i = 0; i < x.size(); i++){
			distance += (x[i] - y[i]) * (x[i] - y[i]);	
		}
		return distance;
	}

	// Check whether CV values are within their respective Voronoi cell in CV space
	bool FiniteTempString::InCell(const CVList& cvs)
	{
		std::vector<double> dists;
		// Record the difference between all cvs and all nodes
		dists.resize(_numnodes);
		for (size_t i = 0; i < _numnodes; i++){
			dists[i] = 0;
			for(size_t j = 0; j < cvs.size(); j++)
				dists[i]+=(cvs[j]->GetValue() - _worldstring[j][i])*(cvs[j]->GetValue() - _worldstring[j][i]);
		}
		
		if(std::min_element(dists.begin(), dists.end()) - dists.begin() != _mpiid){
			return false;
		} else {
			return true;
		}
	}

	// Check whether tolerance criteria has been met
	bool FiniteTempString::TolCheck()
	{
		std::vector<double> toldists;
		double centersize = _centers.size();
		toldists.resize(centersize * _numnodes);

		for(size_t i = 0; i < _numnodes; i++){
			for(size_t j = 0; j < _centers.size(); j++){
					toldists[centersize * i + j] = fabs(_tolcheck[j][i] - _worldstring[j][i]);
			}
		}

		if(*std::max_element(toldists.begin(), toldists.end()) < _tol){
			return true;
		} else {
			return false;
		}

	}

	// Pre-simulation hook
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

		_iterator = 0;

		// Used for reparameterization
	 	_alpha = _mpiid / (_numnodes - 1.0);

		// Initialize running averages
		for(size_t i = 0; i< _centers.size(); ++i){
			_worldstring[i].resize(_numnodes);
			_runavgs[i] = 0;
			_cv_prev[i] = cvs[i]->GetValue();

			// Gathers into _worldstring, where it is cv index followed by node index
			mpi::all_gather(_world, _centers[i], _worldstring[i]);
		}

		if(_restart){
			_currentiter = _restartiter + 1;
			for(size_t i = 0; i < _centers.size(); i++){
				_runavgs[i] = _restartavgs[i];
			}
			_stringout << "#! Restarting FTS from previous run";
		} else {
			_stringout << "#! Beginning new FTS run";
		}

		_stringout << " with " << _numnodes << " nodes starting at iteration " << _currentiter << std::endl;

		if(_tol){
			_stringout << "#! Tolerance criteria is ON: " << _tol << std::endl;
			_tolcheck.resize(_centers.size());
			for(size_t i = 0; i < _centers.size(); i++){
				_tolcheck[i].resize(_numnodes);
			}
		} else {
			_stringout << "#! Tolerance criteria is OFF" << std::endl;
		}

		if(_maxiterator){
			_stringout << "#! Max iterations is ON: " << _maxiterator << std::endl;
		} else {
			_stringout << "#! Max iterations is OFF" << std::endl;
		}

		_stringout << "#! Walker# Iteration#";
		for(size_t i = 0; i < _centers.size(); i++){
			_stringout << " NodeCV" << i << " CurrentValueCV" << i << " RunningAvgCV" << i;
		}
		_stringout << std::endl;

		for(size_t i=0; i<positions.size();i++){
			for(size_t j=0; j<positions[i].size();j++)
				_prev_positions[i][j] = positions[i][j];
		}

		// If initial config. is not in Voronoi cell, apply umbrella potential to move into the cell
		_run_SMD = !InCell(cvs);

	}

	// Post-integration hook.
	void FiniteTempString::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{

		auto& forces = snapshot->GetForces();
		auto& positions = snapshot->GetPositions();

		// Apply umbrella potential in increments of 2000 time steps
		if(_run_SMD && _cv_inside_iterator == 2000){
			if(InCell(cvs)){
				_run_SMD = false;
				for(auto& force : forces){
					for(auto& xyz : force){
						xyz = 0.0;
					}
				}
			} else {
				_cv_inside_iterator = 0;
			}
		}
		if(_run_SMD){
			for(size_t i = 0; i < cvs.size(); i++){
				// Get current cv and gradient
				auto& cv = cvs[i];
				auto& grad = cv->GetGradient();

				// Compute dV/dCV
				auto D = _spring*(cv->GetDifference(_centers[i]));

				// Update forces
				for(size_t j = 0; j < forces.size(); j++){
					for(size_t k = 0; k < forces[j].size(); k++){
						forces[j][k] -= D*grad[j][k];
					}
				}
			}
			
			_cv_inside_iterator++;
		}
		else
		{
			auto inside = InCell(cvs);
			// If new step is within the Voronoi cell, continue with FTS sampling
			// Otherwise, reset system back to previous step  
			if(!inside){
				for(auto& force : forces)
					for(auto& xyz : force)
						xyz = 0.0;

				for(size_t i = 0; i < positions.size(); i++)
					for(size_t j = 0; j < positions[i].size(); j++)
						positions[i][j] = _prev_positions[i][j];	
			}
			else{
				for(size_t i = 0; i < positions.size(); i++)
					for(size_t j = 0; j < positions[i].size(); j++)
						_prev_positions[i][j] = positions[i][j];

				for(size_t i =0; i < cvs.size(); i++)
					_cv_prev[i] = cvs[i]->GetValue();
			}

			// Calculate running averages for each CV at each node 
			for(size_t i = 0; i < _runavgs.size(); i++){
				_runavgs[i] = _runavgs[i] * (_currentiter * _blockiterations + _iterator) + _cv_prev[i];
				_runavgs[i] /= (_currentiter * _blockiterations + _iterator + 1);
			}

			// Update the string, every _blockiterations string method iterations
			if(_iterator > _blockiterations){

				// Write out the string to file
				PrintString(cvs);

				// Update the string and reparameterize 
				StringUpdate();

				_iterator = 0;
				_currentiter++;

				if(_maxiterator && _currentiter > _maxiterator){
					std::cout << "System has reached max string method iterations (" << _maxiterator << ") as specified in the input file(s)." << std::endl;
					std::cout << "Exiting now" << std::endl;
					_world.abort(-1);
				}

				for(size_t i = 0; i < _centers.size(); i++)
					mpi::all_gather(_world, _centers[i], _worldstring[i]);

				if(_tol){
					// check whether max distance moved in any CV is below tolerance criteria _tol
					if(TolCheck()){
						std::cout << "System has converged within tolerance criteria of " << _tol << std::endl;
						std::cout << "Exiting now" <<std::endl;
						_world.abort(-1);
					}
				}

				if(!InCell(cvs)){
					_run_SMD = true;
					_cv_inside_iterator = 0;
				}
			} else {
				_iterator++;
			}
		}
	}

	// Post-simulation hook.
	void FiniteTempString::PostSimulation(Snapshot*, const CVList&)
	{
		_stringout.close();
		std::cout << "Walker " << _mpiid << " has reached max MD steps as specified in the input file(s)." << std::endl;
		std::cout << "Completed iterations: " << _currentiter << std::endl;
		std::cout << "Exiting now" << std::endl;
		_world.abort(-1);

	}

	void FiniteTempString::PrintString(const CVList& CV)
	{
		_stringout.precision(8);
		_stringout << _mpiid << " " << _currentiter << " ";
		for(size_t i = 0; i < _centers.size(); i++){
			_stringout << _centers[i] << " " << CV[i]->GetValue() << " " << _runavgs[i] << " ";
		}
		_stringout << std::endl;
	}

	void FiniteTempString::StringUpdate()
	{
		size_t i;
		size_t centersize = _centers.size();
		double alpha_star;
		std::vector<double> alpha_starv;
		std::vector<double> cvs_new;
		std::vector<double> cvs_newv; 

		std::vector<double> lcv0, ucv0;
		lcv0.resize(centersize);
		ucv0.resize(centersize);

		int sendneighbor, recvneighbor;
		MPI_Status status;

		if(_mpiid == 0){
			sendneighbor = 1;
			recvneighbor = _world.size()-1;
		} 
		else if (_mpiid == (unsigned int)_world.size() - 1){
			sendneighbor = 0;
			recvneighbor = _world.size() - 2;
		} 
		else{
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

		// Update node locations toward running averages:
		for(i = 0; i < cvs_new.size(); i++){
			if(_mpiid == 0 || _mpiid == _numnodes - 1)
				cvs_new[i] = _centers[i] - _tau * (_centers[i] - _runavgs[i]);
			else
				cvs_new[i] = _centers[i] - _tau * (_centers[i] - _runavgs[i]) + 
					(_kappa * _numnodes * _tau * 
					(ucv0[i] + lcv0[i] - 2 * _centers[i]));
		}

		// Reparameterization:
		// alpha_star ranges from 0 to 1 from one end of string to other
		if(_mpiid == 0)
			alpha_star = 0;
		else
			alpha_star = sqrt(sqdist(_centers, lcv0));

		mpi::all_gather(_world, alpha_star, alpha_starv);

		for(i = 1; i < alpha_starv.size(); i++){
			alpha_starv[i] = alpha_starv[i-1] + alpha_starv[i];
		}

		for(i = 1; i < alpha_starv.size(); i++){
			alpha_starv[i] /= alpha_starv[_numnodes - 1];
		}

		tk::spline spl;

		if(_tol){
			for(i = 0; i < _centers.size(); i++){
				mpi::all_gather(_world, _centers[i],_tolcheck[i]);
			}
		}

		for(i = 0; i < centersize; i++){
			mpi::all_gather(_world, cvs_new[i], cvs_newv);
			// spl.setpoints(alpha_starv, vector of all new points **in one particular dimension**);
			spl.set_points(alpha_starv, cvs_newv);
			// Node locations are updated with post-reparameterization values
			_centers[i] = spl(_alpha);
		}

	}
}
