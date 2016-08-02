#include "ABF.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <cassert>


//Later, might be quite worthwhile to adopt the EXTENDED DYNAMICS approach instead
//which works with constraints and multiple CVs easily
//Then, this code can be used to compare results.


//There is a better way by Darve, which uses a time derivative instead of a 2nd order spatial
//
//There are also other ways. The dX/dE terms involve a lot of freedom. (inverse gradients)
//Details supplied below
//
//NOTABLY: Convergence RATE likely depends on the choice of expression for F
//
//The multiple expressions comes from the inverse gradient term, which can be defined in many ways
//The propogation of the force to cartesian coordinates can be done through any vector field
//subject to orthonormality conditions. For details, see above publication and:

//ADD PUBS
//


//I lied. I will use the following definition for a vector field that fits all
//the orthonormality requirements (from above):

//del(E)/|del(E)|^2

//and use that as w in the equation below:

//F = -d/dt(w.p)
//where:
//F is the instantaneous force
//-d/dt is the time derivative
//w is any vector field that satisfies w.grad(E) = 1
//p is momenta

//Equation is from J. Chem. Phys. (2008)
//"Adaptive biasing force method for scalar and vector free energy calculations"
//Darve, Rodriguez-Gomez, Pohorille

//I will use a numerical derivative for time at the same order as the verlet integrator.

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
			if ((_histdetails[i][0] + j*((_histdetails[i][1]-_histdetails[i][0])/_histdetails[i][2]) < cvs[i]->GetValue()) && (_histdetails[i][0] + (j+1)*((_histdetails[i][1]-_histdetails[i][0])/_histdetails[i][2]) > cvs[i]->GetValue()))
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
				{
				indexer = indexer*_histdetails[j][2];
				}
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
		{
		 	_worldout.open(_filename.c_str());
		 }
		
		// Convenience. Number of CVs.
		_dim = cvs.size();

		// Size and initialize _Fold
		_Fold.resize(_dim);
		for(size_t i = 0; i < _Fold.size(); ++i)
			{
			_Fold[i] = 0;
			}
		
		// Size and initialize _Fworld and _Nworld		
		_Fworld.resize(1);
		_Nworld.resize(1);

		for(size_t i = 0; i < _dim; ++i)
			{
			_Fworld.resize(_Fworld.size()*_histdetails[i][2]);
			_Nworld.resize(_Nworld.size()*_histdetails[i][2]);
			}
		_Fworld.resize(_Fworld.size()*_dim);


		// Check if _F, _N is loaded or not
		if (_F.size() == _Fworld.size() && _N.size() == _Nworld.size())
		{		
			_walkerout << "Loaded histogram from user input." << std::endl;
			_walkerout << "Restarting from iteration "<< _iteration <<"." << std::endl;
		}

		else if (_F.size() == 0 && _N.size() == 0)
		{		
			_walkerout << "No user provided initial histogram was found in input. Initializing to zero." << std::endl;
			_F.resize(1);
			_N.resize(1);
			for(size_t i = 0; i < _dim; ++i)
				{
				_F.resize(_F.size()*_histdetails[i][2]);
				_N.resize(_N.size()*_histdetails[i][2]);
				}
			_F.resize(_F.size()*_dim);
		}
		else
		{
			_walkerout << "F or N was specified in input, but failed to load correctly. Please check that dimension and the length are correct. Please refer to the manual to see an example on loading histograms.";
			
			if(_mpiid != 0)
		 		_worldout.open(_filename.c_str());
			_worldout << "F or N was specified in input, but failed to load correctly for at least one walker. Please check that dimension and the length are correct. Refer to walker outputs to find which walker(s) had errors";
			_worldout.close();
				
			std::cerr << "F or N was specified in input, but failed to load correctly for at least one walker. Please check that dimension and the length are correct. Refer to walker outputs to find which walker(s) had errors";
			_walkerout.close();			
			exit(EXIT_FAILURE);
		}
		

		// Janky way to set the size right.
		_biases.resize((cvs[0]->GetGradient()).size());
		for(size_t i = 0; i < _biases.size(); ++i)
			{
			_biases[i].resize(3);
			for(size_t j = 0; j < _biases[i].size(); ++j)
				{
				_biases[i][j] = 0;
				}
			}
		
		//One wdotpold for each CV
		_wdotpold.resize(2*cvs.size());
		for(size_t i = 0; i < _wdotpold.size(); ++i)
			{
			_wdotpold[i] = 0.0;
			}

		// Print out the coordinate system for _F and _N
		if(_mpiid == 0)
			{
			_worldout << "Histogram grids are set up. The coordinates for Generalized Force Vector Field and N, in order of CVs in columns are:";

			int modulo = 1;
			int index;
			for(size_t i = 0; i < _N.size(); ++i)
				{
				index = i;
				_worldout << std::endl;
				for(size_t j=0 ; j < _histdetails.size(); ++j)
					{
					modulo = 1;
					for(size_t k=j+1 ; k <_histdetails.size(); ++k)
						{
						modulo = modulo * _histdetails[k][2];
						}
					_worldout << (floor(index/modulo)+0.5)*((_histdetails[j][1]-_histdetails[j][0])/_histdetails[j][2]) + _histdetails[j][0] << " ";
					index = index % modulo;
					}
				}
			}
		
		_walkerout << std::endl;
		_walkerout << std::endl;
		_walkerout << "Start of simulation. Printing out 1) Dimensionality of F 2) Total Length of N 3) Total Length of F, 4) CV details 5) Initial histogrammed F estimate." << std::endl;
		_walkerout << _dim << std::endl;
		_walkerout << _N.size() <<std::endl;
		_walkerout << _F.size() <<std::endl;
		_walkerout << "CV bounds and nr of bins are: ";
		for(size_t i = 0; i < _histdetails.size(); ++i)
			{
			_walkerout << "[" << _histdetails[i][0] << "," << _histdetails[i][1] << "," << _histdetails[i][2] << "] ";
			}
		_walkerout << std::endl;

		// Print out the initial estimate. User can verify whether their provided estimate was successfully loaded.
		_walkerout << "The initial F and N for this walker is printed below. If loaded from an input, please check for accuracy." << std::endl;
		for(size_t i = 0; i < _N.size(); ++i)
				{
				for(size_t j = 0; j < _dim; ++j)
					{
					_walkerout << _F[_dim*i+j] << " ";
					}
				_walkerout << _N[i] << std::endl;
				}
			_walkerout << std::endl;
		
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
		// Feedback on Iteration Number. Helps with subsequent prints and comparing with log files.
		if(_iteration % _printdetails[0] == 0)
			_walkerout << "Iteration: " <<_iteration;		
		
		// One genforce per CV.		
		std::vector<double> genforce(cvs.size());

		// Gather information. Can clean this up later.
		auto& vels = snapshot->GetVelocities();
		auto& mass = snapshot->GetMasses();
		auto& forces = snapshot->GetForces();
		auto& positions = snapshot->GetPositions();
		auto& ids = snapshot->GetAtomIDs();
		auto& types = snapshot->GetAtomTypes();

		//! Coord holds where we are in CV space in current timestep.
		int coord = histCoords(cvs);

		//! Variable to hold w dot p scalar for each CV dimension. Since memory of this is not required, will reset to 0 for every CV.
		double wdotp = 0;

		
		std::vector<double> dellensq(cvs.size()); /**< Vector to keep track of the norm of projector for normalization. Dellensq will start at 0 at the beginning of the loop for each CV, and the values of the gradient squared will be added through the loop. At the end of the loop, it will be equal to |grad(CV)|^2, and can be used to normalize the vector field. Is also used for Gram-Schmidt Orthogonalization, if enabled. */

		//! Vector<vector<double>> to hold the projector vector field.
		/*! 
		 * Projector holds w, the projection of the force from generalized coordinates to cartesian. 
		 * To save a dimension, it holds x, y and z sequentially, and is thus 3xNrAtoms long for each CV.
		 * If Gram-Schmidt Orthogonalization is enabled, the first projector is grad(CV)/|grad(CV)|^2 and the rest are derived from the algorithm.
		 * If GSO is disabled, each projector is simply grad(CV)/|grad(CV)|^2.
		 */
		std::vector<std::vector<double>> projector(cvs.size(), std::vector<double>(vels.size()*vels[0].size()));

		//Bound check on CV, implemented in histCoords function. If in bounds, bias and collect values for histogram.
		if (coord != -1)
			{
			// Add to histogram hits.
			_N[coord] += 1;
		
			// Loop to calculate generalized force for each CV
			// Outmost loop through CVs
			for(size_t i = 0; i < cvs.size(); ++i)
				{
				//Feedback on CVs.
				if(_iteration % _printdetails[0] == 0 && _printdetails[1] > 0)
					_walkerout << " CV[" << i << "]: " << cvs[i]->GetValue();
				
				// Get the gradient of current CV.
				auto& grad = cvs[i]->GetGradient();
				
				// Reset w dot p
				wdotp = 0;
			
				// Reset dellensq
				dellensq[i] = 0;

				if(_Orthogonalization) // Begin Gram-Schmidt Orthogonalization here. For all but the first CV, will remove the projection of every other CVs vector field from the current CV projector.
					{
					if(_iteration % _printdetails[0] == 0 && _printdetails[2] > 0 && i > 0)
						_walkerout << " Orthogonalization Corrector vector is: ";
					for(size_t j = 0; j < i; ++j)
						{
						for(size_t k = 0; k < grad.size(); ++k)
							{
							for(size_t l = 0; l < grad[k].size(); ++l)
								{
								projector[i][3*k+l] -= grad[k][l]*projector[j][3*k+l]/dellensq[j];
								if(_iteration % _printdetails[0] == 0 && _printdetails[2] > 0 && grad[k][l] != 0)
									_walkerout <<  projector[i][3*k+l] << " ";
								}
							}
						}
					} // End Gram-Schmidt Orthogonalization. Projector should be zero for each entry if a) This is disabled, or b) if CVs are orthogonal.
				
				// Start calculation of w dot p
				for(size_t j = 0; j < grad.size(); ++j)
					{
					for(size_t k = 0; k <grad[j].size(); ++k)
						{
						projector[i][3*j+k] += grad[j][k];			
						wdotp += (projector[i][3*j+k])*vels[j][k]*mass[j];
						dellensq[i] += pow(projector[i][3*j+k],2);						
						}		
					}
			
				// Normalize to ensure JW = I
				wdotp = wdotp/dellensq[i];
				
				// Feedback on dellensq.
				if(_iteration % _printdetails[0] == 0 && _printdetails[3] > 0)
					_walkerout << " Normalization factor is: " << dellensq[i];
	
				// Calculate instantaneous generalized force. This is a second order forward numerical time derivative. +_Fold[i] removes the bias.
				genforce[i] = _unitconv*(3*wdotp-4*_wdotpold[i]+_wdotpold[i+cvs.size()])/(2*_timestep) + _Fold[i];

				// Add to running total in current bin.
				_F[_dim*coord+i] += genforce[i];

				// Sync _F and _N across all walkers.
				MPI_Allreduce(&_N[0], &_Nworld[0], _N.size(), MPI_INT, MPI_SUM, _world);
				MPI_Allreduce(&_F[0], &_Fworld[0], _F.size(), MPI_DOUBLE, MPI_SUM, _world);
		
				// Keep track of bias to remove it in next timestep.
				_Fold[i] = _Fworld[_dim*coord+i]/std::max(_min,_Nworld[coord]);

				// Record t=i-1 and t=i-2 values of w dot p for time derivative.
				_wdotpold[i+cvs.size()] = _wdotpold[i];
				_wdotpold[i] = wdotp;
				}		
			}
		else //Bound check on CV, implemented in histCoords function. If outside bounds, do not bias and do not collect histogram, but keep track of wdotp.
			{
			// Feedback on a CV being out of bounds.
			if(_iteration % _printdetails[0] == 0 && _printdetails[1] > 0)
				_walkerout << " At least one CV is out of bounds!";

			for(size_t i = 0; i < cvs.size(); ++i)
				{
				// Feedback on CV.
				if(_iteration % _printdetails[0] == 0 && _printdetails[1] > 0)
					_walkerout << " CV[" << i << "]: " << cvs[i]->GetValue();

				wdotp = 0;
				dellensq[i] = 0;		
				auto& grad = cvs[i]->GetGradient();
				if(_Orthogonalization) // Begin Gram-Schmidt Orthogonalization here. For all but the first CV, will remove the projection of every other CVs vector field from the current CV projector.
					{
					if(_iteration % _printdetails[0] == 0 && _printdetails[2] > 0 && i > 0)
						_walkerout << " Orthogonalization Corrector vector is: ";
					for(size_t j = 0; j < i; ++j)
						{
						for(size_t k = 0; k < grad.size(); ++k)
							{
							for(size_t l = 0; l < grad[k].size(); ++l)
								{
								projector[i][3*k+l] -= grad[k][l]*projector[j][3*k+l]/dellensq[j];
								if(_iteration % _printdetails[0] == 0 && _printdetails[2] > 0 && grad[k][l] != 0)
									_walkerout <<  projector[i][3*k+l] << " ";
								}
							}
						}
					} // End Gram-Schmidt Orthogonalization. Projector should be zero for each entry if a) This is disabled, or b) if CVs are orthogonal.
				
				// Start calculation of w dot p
				for(size_t j = 0; j < grad.size(); ++j)
					{
					for(size_t k = 0; k <grad[j].size(); ++k)
						{
						projector[i][3*j+k] += grad[j][k];			
						wdotp += (projector[i][3*j+k])*vels[j][k]*mass[j];
						dellensq[i] += pow(projector[i][3*j+k],2);						
						}		
					}

				
				// Normalize to ensure JW = I
				wdotp = wdotp/dellensq[i];
				genforce[i] = _unitconv*(3*wdotp-4*_wdotpold[i]+_wdotpold[i+cvs.size()])/(2*_timestep) + _Fold[i];
				_Fold[i] = 0;
				_wdotpold[i+cvs.size()] = _wdotpold[i];
				_wdotpold[i] = wdotp;
			}
		}
		
		// Print out current F estimate over CV space - for backup and troubleshooting
		if(_iteration % _FBackupInterv == 0 && _mpiid == 0)
			{
			_worldout << std::endl;
			_worldout << "Iteration: " << _iteration << std::endl;			
			_worldout << "Printing out the current Adaptive Biasing Vector Field." << std::endl;
			_worldout << "First (Nr of CVs) columns are the coordinates, the next (Nr of CVs) columns are components of the Adaptive Force vector at that point." << std::endl;
			_worldout << "The columns are " << _N.size() << " long, mapping out a surface of ";
			for(size_t i = 0; i < _histdetails.size()-1; ++i)
				{
				_worldout << _histdetails[i][2] << " by ";
				}
			_worldout << _histdetails[_histdetails.size()-1][2] <<" points in " << _histdetails.size() << " dimensions." <<std::endl;


			int modulo = 1;
			int index;
			for(size_t i = 0; i < _Nworld.size(); ++i)
				{
				index = i;
				_worldout << std::endl;
				for(size_t j=0 ; j < _histdetails.size(); ++j)
					{
					modulo = 1;
					for(size_t k=j+1 ; k <_histdetails.size(); ++k)
						{
						modulo = modulo * _histdetails[k][2];
						}
					_worldout << (floor(index/modulo)+0.5)*((_histdetails[j][1]-_histdetails[j][0])/_histdetails[j][2]) + _histdetails[j][0] << " ";
					index = index % modulo;
					}
				for(size_t j=0 ; j < _dim; ++j)
					{
					_worldout << _Fworld[_dim*i+j]/std::max(_Nworld[i],_min) << " ";
					}
				}
			_worldout << std::endl;
			_worldout << std::endl;
			}
		
		// Calculate the bias from averaged F at current CV coordinates
		CalcBiasForce(cvs,genforce, snapshot);
		

		// Feedback on biasing.
		if(_iteration % _printdetails[0] == 0 && _printdetails[8] > 0)
			_walkerout << " Biases to forces are: ";

		// Update the forces in snapshot by adding in the force bias from each
		// CV to each atom based on the gradient of the CV.
		for (size_t j = 0; j < forces.size(); ++j)
			{
			for(size_t k = 0; k < forces[j].size(); ++k)
				{
				if(_iteration % _printdetails[0] == 0 && _printdetails[8] > 0 && _biases[j][k] != 0)
					_walkerout << _biases[j][k] << " ";
				forces[j][k] += _biases[j][k];					
				}
			}
		
		// Feedback on locations.
		if(_iteration % _printdetails[0] == 0 && _printdetails[6] > 0)
			{
			_walkerout << " Locations are: ";
			for(size_t i = 0; i < cvs.size(); ++i)
				{
				auto& grad = cvs[i]->GetGradient();
				auto& loc = snapshot->GetPositions();
				for (size_t j = 0; j < forces.size(); ++j)
					{
					for(size_t k = 0; k < forces[j].size(); ++k)
						{
						if (grad[j][k] != 0)
							{
							_walkerout << loc[j][k] << " ";
							}
						}
					}
				}
			}
			
		
	// Each output gets its own line
	if(_iteration % _printdetails[0] == 0)
		_walkerout << std::endl;	
	}

	// Post-simulation hook.
	void ABF::PostSimulation(Snapshot*, const CVList&)
	{
		if(_mpiid == 0)
			{
			_worldout << "End of simulation. The coordinates for Generalized Force Vector Field and N, in order of CVs in columns are:";
			_worldout << std::endl;
			_worldout << "Iteration: " << _iteration << std::endl;			
			_worldout << "Printing out the current Adaptive Biasing Vector Field." << std::endl;
			_worldout << "First (Nr of CVs) columns are the coordinates, the next (Nr of CVs) columns are components of the Adaptive Force vector at that point." << std::endl;
			_worldout << "The columns are " << _N.size() << " long, mapping out a surface of ";
			for(size_t i = 0; i < _histdetails.size()-1; ++i)
				{
				_worldout << _histdetails[i][2] << " by ";
				}
			_worldout << _histdetails[_histdetails.size()-1][2] <<" points in " << _histdetails.size() << " dimensions." <<std::endl;


			int modulo = 1;
			int index;
			for(size_t i = 0; i < _Nworld.size(); ++i)
				{
				index = i;
				_worldout << std::endl;
				for(size_t j=0 ; j < _histdetails.size(); ++j)
					{
					modulo = 1;
					for(size_t k=j+1 ; k <_histdetails.size(); ++k)
						{
						modulo = modulo * _histdetails[k][2];
						}
					_worldout << (floor(index/modulo)+0.5)*((_histdetails[j][1]-_histdetails[j][0])/_histdetails[j][2]) + _histdetails[j][0] << " ";
					index = index % modulo;
					}
				for(size_t j=0 ; j < _histdetails.size(); ++j)
					{
					_worldout << _Fworld[i+j*_Nworld.size()]/_Nworld[i] << " ";
					}
				}
			_worldout.close();
			}				
		_walkerout.close();
	}

	void ABF::CalcBiasForce(const CVList& cvs, const std::vector<double>& genforce, const Snapshot* snapshot)
	{	
		// Reset the bias.
		for(size_t i = 0; i < _biases.size(); ++i)
			for(size_t j = 0; j < _biases[i].size(); ++j)
				_biases[i][j] = 0;
			
		// Figure out where we are in CV space in current timestep
		int coord = histCoords(cvs);

		// Calculate the bias if within CV bounds.
		if (coord != -1)
			{
			for(size_t i = 0; i < cvs.size(); ++i)
				{
				auto& grad = cvs[i]->GetGradient();
				//_walkerout << " Gradient[" << i << "]: ";
				for(size_t j = 0; j < _biases.size(); ++j)
					{
					for(size_t k = 0; k < _biases[j].size(); ++k)
						{
						if (grad[j][k] != 0)
							{
							//_walkerout << grad[j][k] << " ";
							_biases[j][k] += -(_Fworld[_dim*coord+i]/std::max(_Nworld[coord],_min))*grad[j][k];
							}
						}
					}			
				}
			}

		else  // Bias according to a harmonic if outside CV bounds, and if this feature is enabled.
		{
			for(size_t i = 0; i < cvs.size(); ++i)
			{
				if(cvs[i]->GetValue() < _restraint[i][0] && _restraint[i][2] > 0)
				{
					if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[7] > 0)
						_walkerout << " Minimum Harmonic Restraint Active on CV: " << i << "; Restraint force is: ";
					auto& grad = cvs[i]->GetGradient();
					for(size_t j = 0; j < _biases.size(); ++j)
					{
						for(size_t k = 0; k < _biases[j].size(); ++k)
						{
							if (grad[j][k] != 0)
							{
								_biases[j][k] += grad[j][k]*(_restraint[i][0] - cvs[i]->GetValue())*_restraint[i][2];
								if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[7] > 0)
									_walkerout << _biases[j][k] << " ";
							}
						}
					}

				}

				else if(cvs[i]->GetValue() > _restraint[i][1] && _restraint[i][2] > 0)
				{
					if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[7] > 0)
						_walkerout << " Maximum Harmonic Restraint Active on CV: " << i << "; Restraint force is: ";

					auto& grad = cvs[i]->GetGradient();
					for(size_t j = 0; j < _biases.size(); ++j)
					{
						for(size_t k = 0; k < _biases[j].size(); ++k)
						{
							if (grad[j][k] != 0)
							{
								_biases[j][k] += grad[j][k]*(_restraint[i][1] - cvs[i]->GetValue())*_restraint[i][2];
								if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[7] > 0)
									_walkerout << _biases[j][k] << " ";
							}
						}
					}	
				}
			}
		}
	}
}



































