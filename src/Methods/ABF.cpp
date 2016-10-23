/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Emre Sevgen <sesevgen@uchicago.edu>
 *                Hythem Sidky <hsidky@nd.edu>
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
#include "ABF.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <cassert>


// "Adaptive biasing force method for scalar and vector free energy calculations"
// Darve, Rodriguez-Gomez, Pohorille
// J. Chem. Phys. (2008)

namespace mpi = boost::mpi;
namespace SSAGES
{
	// Function to return histogram coordinates, given CV values. Returns -1 if any CV value is outside bounds.
	// Regardless of number of CVs, the histogram is a one dimensional object.
	// Therefore, this function works in two parts:
	// 1) Mapping from CV coordinates to an N dimensional fictitious object where N is the number of CVs.
	// -> This fictitious object is (Number of CV bins) large in each dimension.
	// -> Example: If you have two CVs, X and Y, and you bin X into three bins and Y into five bins, the fictitious object is a 2 dimensional 3x5 matrix. 
	
	// 2) Mapping from N dimensions to one dimension.
	// -> Simply constructed by writing out all members of the fictitious object in a line. This is a vector of length equal to the product of number of bins in each dimension.
	// -> For the above example, this is 3.5 = 15 members. First 5 members are Y = [1 2 3 4 5] bins with X = 1. Next 5 members are Y = [1 2 3 4 5] with X = 2 ....

	//! Remove a column from the histogram.
	/*!
	 * \param matrix Matrix storing the histogram.
	 * \param colToRemove Index of the column to be removed.
	 */
	void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
	{
	    unsigned int numRows = matrix.rows();
	    unsigned int numCols = matrix.cols()-1;

	    if( colToRemove < numCols )
		matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

	    matrix.conservativeResize(numRows,numCols);
	}

	int ABF::histCoords(const CVList& cvs)
	{
		// Histogram details: This is a 2 Dimensional object set up as the following:
		// Histdetails is a vector of three vectors, each of those three vectors are (Number of CVs) long.
		// First of these vectors hold the lower bound for the CVs in order.
		// Second vector holds the upper bound for the CVs in order.
		// Third vector holds number of histogram bins for the CVs in order.
	
		// The first mapping starts here.
		for(size_t i = 0; i < _histdetails.size(); ++i)
		{
			// Check if CV is in bounds for each CV.
			if ((cvs[i]->GetValue()<_histdetails[i][0]) || (cvs[i]->GetValue()>_histdetails[i][1]))
				return -1;
		}

		// This vector holds the CV coordinates in the fictitious object.
		std::vector<int> coords(cvs.size());
		
		// Loop over each CV dimension.
		for(size_t i = 0; i < _histdetails.size(); ++i)
		{
			// Loop over all bins in each CV dimensions
			for(int j = 0; j < _histdetails[i][2] ; ++j)
			{
				if((_histdetails[i][0] + j*((_histdetails[i][1]-_histdetails[i][0])/_histdetails[i][2]) < cvs[i]->GetValue()) && 
				   (_histdetails[i][0] + (j+1)*((_histdetails[i][1]-_histdetails[i][0])/_histdetails[i][2]) > cvs[i]->GetValue()))
				{
					coords[i] = j;
				}
			}
		} //First mapping ends here. 

		// Second mapping starts here.
		int finalcoord = 0;
		int indexer = 1;
		// Loop over each CV dimension.
		for(size_t i = 0; i < _histdetails.size(); ++i)
		{
			indexer = 1;
			for(size_t j = i+1; j < _histdetails.size(); ++j)
				indexer = indexer*_histdetails[j][2];

			finalcoord += coords[i]*indexer;
		}

		return finalcoord;			
	}
	
	// Pre-simulation hook.
	void ABF::PreSimulation(Snapshot* snapshot, const CVList& cvs)
	{
		_mpiid = snapshot->GetWalkerID();
		char file[1024];
		sprintf(file, "node-%04d.log",_mpiid);
	 	_walkerout.open(file);

		if(_mpiid == 0)
		 	_worldout.open(_filename.c_str());

		// Convenience. Number of CVs.
		_ncv = cvs.size();

		// Size and initialize _Fold
		_Fold.setZero(_ncv);
		
		// Size and initialize _Fworld and _Nworld		
		auto nel = 1;
		for(auto i = 0; i < _ncv; ++i)
			nel *= _histdetails[i][2];

		_Fworld.setZero(nel*_ncv);
		_Nworld.resize(nel, 0);

		// If F or N are empty, size appropriately. 
		if(_F.size() == 0) _F.setZero(nel*_ncv);
		if(_N.size() == 0) _N.resize(nel, 0);

		// Initialize biases.
		_biases.resize(snapshot->GetPositions().size(), Vector3{0, 0, 0});
		
		// Initialize w \dot p's for finite difference. 
		_wdotp1.setZero(_ncv);
		_wdotp2.setZero(_ncv);
		
	}

	//! Post-integration hook.
	/*!
	 * Post-integration hook is where processes carried out every timestep happen.
	 * First, information from current snapshot is retrieved and stored in variables as necessary.
	 * Then, coordinates in CV space are determined.
	 * Then, for each CV, the time derivative of w.p is calculated.
	 * Then, information is printed out if its enabled.
	 * Finally, bias is applied from current estimate of generalized force.
	 */
	void ABF::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		++_iteration;

		// Gather information.
		auto& vels = snapshot->GetVelocities();
		auto& mass = snapshot->GetMasses();
		auto& forces = snapshot->GetForces();
		auto n = snapshot->GetNumAtoms();

		//! Coord holds where we are in CV space in current timestep.
		int coord = histCoords(cvs);

		//! Eigen::MatrixXd to hold the CV gradient.
		Eigen::MatrixXd J(_ncv, 3*n);

		// Fill J. Each column represents grad(CV) with flattened Cartesian elements. 
		for(int i = 0; i < _ncv; ++i)
		{
			auto& grad = cvs[i]->GetGradient();
			for(size_t j = 0; j < n; ++j)
				J.block<1, 3>(i,3*j) = grad[j];
		}
		
		//* Calculate W using Darve's approach (http://mc.stanford.edu/cgi-bin/images/0/06/Darve_2008.pdf).
		Eigen::MatrixXd Jmass = J.transpose();
		
		for(size_t i = 0; i < forces.size(); ++i)
			Jmass.block(3*i, 0, 3, _ncv) = Jmass.block(3*i, 0, 3, _ncv)/mass[i];
							
		Eigen::MatrixXd Minv = J*Jmass;
		MPI_Allreduce(MPI_IN_PLACE, Minv.data(), Minv.size(), MPI_DOUBLE, MPI_SUM, _comm);

		Eigen::MatrixXd Wt = Minv.inverse()*Jmass.transpose();	

		// Fill momenta.
		Eigen::VectorXd momenta(3*vels.size());
		for(size_t i = 0; i < vels.size(); ++i)
			momenta.segment<3>(3*i) = mass[i]*vels[i];

		// Compute dot(w,p).
		Eigen::VectorXd wdotp = Wt*momenta;

		// Reduce dot product across processors.
		MPI_Allreduce(MPI_IN_PLACE, wdotp.data(), wdotp.size(), MPI_DOUBLE, MPI_SUM, _comm);

		// Compute d(wdotp)/dt second order backwards finite difference. 
		// Adding old force removes bias. 
		Eigen::VectorXd dwdotpdt = _unitconv*(1.5*wdotp - 2.0*_wdotp1 + 0.5*_wdotp2)/_timestep + _Fold;

		// If we are in bounds, sum force into running total.
		if(coord != -1)
		{
			_F.segment(_ncv*coord, _ncv) += dwdotpdt;
			++_N[coord];
		}

		// Reduce data across processors.
		MPI_Allreduce(_F.data(), _Fworld.data(), _F.size(), MPI_DOUBLE, MPI_SUM, _world);
		MPI_Allreduce(_N.data(), _Nworld.data(), _N.size(), MPI_INT, MPI_SUM, _world);

		// If we are in bounds, store the old summed force.
		if(coord != -1)
			_Fold = _Fworld.segment(_ncv*coord, _ncv)/std::max(_min, _Nworld[coord]);
	
		// Update finite difference time derivatives.
		_wdotp2 = _wdotp1;
		_wdotp1 = wdotp;		

		// Write out data to file.
		if(_iteration % _FBackupInterv == 0)
			WriteData();
		
		// Calculate the bias from averaged F at current CV coordinates
		CalcBiasForce(snapshot, cvs, coord);
		
		// Update the forces in snapshot by adding in the force bias from each
		// CV to each atom based on the gradient of the CV.
		for (size_t j = 0; j < forces.size(); ++j)
			forces[j] += _biases[j];	
	}

	// Post-simulation hook.
	void ABF::PostSimulation(Snapshot*, const CVList&)
	{
		WriteData();
		_worldout.close();		
		_walkerout.close();
	}

	void ABF::CalcBiasForce(const Snapshot* snapshot, const CVList& cvs, int coord)
	{	
		// Reset the bias.
		for(auto& b : _biases)
			b.setZero();
		
		// Compute bias if within bounds
		if(coord != -1)
		{
			for(int i = 0; i < _ncv; ++i)
			{
				auto& grad = cvs[i]->GetGradient();
				for(size_t j = 0; j < _biases.size(); ++j)
					_biases[j] -= _Fworld[_ncv*coord+i]*grad[j]/std::max(_min, _Nworld[coord]);
			}
		}
		else
		{
			for(int i = 0; i < _ncv; ++i)
			{
				auto k = 0.;
				auto x0 = 0.;

				if(cvs[i]->GetValue() < _restraint[i][0] && _restraint[i][2] > 0)
				{
					k = _restraint[i][2];
					x0 = _restraint[i][0];
				}
				else if (cvs[i]->GetValue() > _restraint[i][1] && _restraint[i][2] > 0)
				{
					k = _restraint[i][2];
					x0 = _restraint[i][1];
				}

				auto& grad = cvs[i]->GetGradient();
				for(size_t j = 0; j < _biases.size(); ++j)
					_biases[j] -= (cvs[i]->GetValue() - x0)*k*grad[j];
			}	
		}
	}

	void ABF::WriteData()
	{		
		if(_mpiid != 0)
			return;

		_worldout << std::endl;
		_worldout << "Iteration: " << _iteration << std::endl;			
		_worldout << "Printing out the current Adaptive Biasing Vector Field." << std::endl;
		_worldout << "First (Nr of CVs) columns are the coordinates, the next (Nr of CVs) columns are components of the Adaptive Force vector at that point." << std::endl;
		_worldout << "The columns are " << _N.size() << " long, mapping out a surface of ";
		
		for(size_t i = 0; i < _histdetails.size()-1; ++i)
			_worldout << _histdetails[i][2] << " by ";
		
		_worldout << _histdetails[_histdetails.size()-1][2] 
		          << " points in " << _histdetails.size() 
		          << " dimensions." <<std::endl;

		int modulo = 1;
		int index = 0;
		for(size_t i = 0; i < _Nworld.size(); ++i)
		{
			index = i;
			_worldout << std::endl;
			for(size_t j=0 ; j < _histdetails.size(); ++j)
			{
				modulo = 1;
				for(size_t k=j+1 ; k <_histdetails.size(); ++k)
					modulo = modulo * _histdetails[k][2];

				_worldout << (floor(index/modulo)+0.5)*((_histdetails[j][1]-_histdetails[j][0])/_histdetails[j][2]) + _histdetails[j][0] << " ";
				index = index % modulo;
			}

			for(int j = 0; j < _ncv; ++j)
				_worldout << _Fworld[_ncv*i+j]/std::max(_Nworld[i],_min) << " ";
		}

		_worldout << std::endl;
		_worldout << std::endl;
	}
}



































