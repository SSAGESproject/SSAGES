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
	 	_stringout.open(file);
		// Convenience. Number of CVs.
		_dim = cvs.size();

		// Size and initialize _Fold
		_Fold.resize(_dim);
		for(size_t i = 0; i < _Fold.size(); ++i)
			_Fold[i] = 0;

		// Size and initialize _F and _N		
		_F.resize(1);
		_N.resize(1);
		
		_Fworld.resize(1);
		_Nworld.resize(1);
		for(size_t i = 0; i < _dim; ++i)
			{
			_F.resize(_F.size()*_histdetails[i][2]);
			_N.resize(_N.size()*_histdetails[i][2]);

			_Fworld.resize(_Fworld.size()*_histdetails[i][2]);
			_Nworld.resize(_Nworld.size()*_histdetails[i][2]);
			}

		_F.resize(_F.size()*_dim);
		_Fworld.resize(_Fworld.size()*_dim);
		


		// readF == " " means do not load a histogram, start from scratch.
		if (_readF == " ")
			{
			_stringout << "No input file specified. Starting from an empty histogram." << std::endl;
			for(size_t i = 0; i < _N.size(); ++i)
				{
				for(size_t j = 0; j < _dim; ++j)
					{
					_F[_dim*i+j] = 0;
					}
				_N[i] = 0;
				}
			}
		else 
			{
			_stringout << "There is an input file specified. Input file is: " << _readF << std::endl;
			std::fstream input(_readF, std::ios_base::in);
			if(input.is_open())
				{
				_stringout << "Loading F estimate from file." <<std::endl;
				unsigned int dimFile;
				input >> dimFile;
				unsigned int len;
				input >> len;
				_stringout << "File dim: "<< dimFile << "; File length: " << len << std::endl;
		
				for(size_t i=0; i < dimFile; ++i)
					input.ignore(256,']');
						

				// Check if fed histogram is the correct dimension and number of bins.
				if (dimFile == _dim && len == _N.size())
					{
					int tempN;
					double tempF;
					for(size_t i = 0; i < _N.size(); ++i)
						{
						for(size_t j = 0; j < _dim; ++j)
							{
							input >> tempF;
							_F[_dim*i+j] = tempF;
							}						
						input >> tempN;
						_N[i] = tempN;
						}
					}
				else
					{
					std::cerr << "An input file was specified, but failed to load correctly. Please check that dimension and the length are correct, and the header is correct. Please refer to the manual to see an example on loading histograms.";
					exit(EXIT_FAILURE);
					//_stringout << "Dimension or Number of Bins mismatch! Initializing Empty Histogram Instead!" << std::endl;
					//for(size_t i = 0; i < _F.size(); ++i)
						//{
						//_F[i].resize(dim);
						//for(size_t j = 0; j < dim; ++j)
						//	_F[i][j] = 0;
						//_N[i] = 0;
						//}
					}
				}
			else
				{
				std::cerr << "An input file was specified, but was not found.";
				exit(EXIT_FAILURE);
				}
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
			_wdotpold[i] = 0.0;
		
		// Print out the coordinate system for _F and _N
		_stringout << "Histogram grids are set up. The coordinates for Generalized Force Vector Field and N, in order of CVs in columns are:";

		int modulo = 1;
		int index;
		for(size_t i = 0; i < _N.size(); ++i)
			{
			index = i;
			_stringout << std::endl;
			for(size_t j=0 ; j < _histdetails.size(); ++j)
				{
				modulo = 1;
				for(size_t k=j+1 ; k <_histdetails.size(); ++k)
					{
					modulo = modulo * _histdetails[k][2];
					}
				_stringout << (floor(index/modulo)+0.5)*((_histdetails[j][1]-_histdetails[j][0])/_histdetails[j][2]) + _histdetails[j][0] << " ";
				index = index % modulo;
				}
			}
		
		_stringout << std::endl;
		_stringout << std::endl;
		_stringout << "Start of simulation. Printing out 1) Dimensionality of F 2) Total Length of N 3) Total Length of F, 4) CV details 5) Initial histogrammed F estimate." << std::endl;
		_stringout << _dim << std::endl;
		_stringout << _N.size() <<std::endl;
		_stringout << _F.size() <<std::endl;
		_stringout << "CV bounds and nr of bins are: ";
		for(size_t i = 0; i < _histdetails.size(); ++i)
			{
			_stringout << "[" << _histdetails[i][0] << "," << _histdetails[i][1] << "," << _histdetails[i][2] << "] ";
			}
		_stringout << std::endl;

		// Print out the initial estimate. User can verify whether their provided estimate was successfully loaded.
		for(size_t i = 0; i < _N.size(); ++i)
				{
				for(size_t j = 0; j < _dim; ++j)
					{
					_stringout << _F[_dim*i+j] << " ";
					}
				_stringout << _N[i] << std::endl;
				}
			_stringout << std::endl;
		
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
	
		// Feedback on Iteration Number. Helps with subsequent prints and comparing with log files.
		if(snapshot->GetIteration() % _printdetails[0] == 0)
			_stringout << "Iteration: " <<snapshot->GetIteration();		
		
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
				if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[1] > 0)
					_stringout << " CV[" << i << "]: " << cvs[i]->GetValue();
				
				// Get the gradient of current CV.
				auto& grad = cvs[i]->GetGradient();
				
				// Reset w dot p
				wdotp = 0;
			
				// Reset dellensq
				dellensq[i] = 0;

				if(_Orthogonalization) // Begin Gram-Schmidt Orthogonalization here. For all but the first CV, will remove the projection of every other CVs vector field from the current CV projector.
					{
					if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[2] > 0 && i > 0)
						_stringout << " Orthogonalization Corrector vector is: ";
					for(size_t j = 0; j < i; ++j)
						{
						for(size_t k = 0; k < grad.size(); ++k)
							{
							for(size_t l = 0; l < grad[k].size(); ++l)
								{
								projector[i][3*k+l] -= grad[k][l]*projector[j][3*k+l]/dellensq[j];
								if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[2] > 0 && grad[k][l] != 0)
									_stringout <<  projector[i][3*k+l] << " ";
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
				if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[3] > 0)
					_stringout << " Normalization factor is: " << dellensq[i];
	
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
			if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[1] > 0)
				_stringout << " At least one CV is out of bounds!";

			for(size_t i = 0; i < cvs.size(); ++i)
				{
				// Feedback on CV.
				if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[1] > 0)
					_stringout << " CV[" << i << "]: " << cvs[i]->GetValue();

				wdotp = 0;
				dellensq[i] = 0;		
				auto& grad = cvs[i]->GetGradient();
				if(_Orthogonalization) // Begin Gram-Schmidt Orthogonalization here. For all but the first CV, will remove the projection of every other CVs vector field from the current CV projector.
					{
					if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[2] > 0 && i > 0)
						_stringout << " Orthogonalization Corrector vector is: ";
					for(size_t j = 0; j < i; ++j)
						{
						for(size_t k = 0; k < grad.size(); ++k)
							{
							for(size_t l = 0; l < grad[k].size(); ++l)
								{
								projector[i][3*k+l] -= grad[k][l]*projector[j][3*k+l]/dellensq[j];
								if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[2] > 0 && grad[k][l] != 0)
									_stringout <<  projector[i][3*k+l] << " ";
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
		if(snapshot->GetIteration() % _FBackupInterv == 0)
		{
		_stringout << std::endl;
		_stringout << "Printing out 1) Dimensionality of F 2) Total Length of N 3) Total Length of F, 4) CV details 5) Initial histogrammed F estimate." << std::endl;
		_stringout << _dim << std::endl;
		_stringout << _Nworld.size() <<std::endl;
		_stringout << _Fworld.size() <<std::endl;
		_stringout << "CV bounds and nr of bins are: ";
		for(size_t i = 0; i < _histdetails.size(); ++i)
			{
			_stringout << "[" << _histdetails[i][0] << "," << _histdetails[i][1] << "," << _histdetails[i][2] << "] ";
			}
		_stringout << std::endl;
		for(size_t i = 0; i < _Nworld.size(); ++i)
				{
				for(size_t j = 0; j < _dim; ++j)
					{
					_stringout << _Fworld[_dim*i+j] << " ";
					}
				_stringout << _Nworld[i] << std::endl;
				}
			_stringout << std::endl;
		}
		
		// Calculate the bias from averaged F at current CV coordinates
		CalcBiasForce(cvs,genforce, snapshot);
		

		// Feedback on biasing.
		if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[8] > 0)
			_stringout << " Biases to forces are: ";


		// Update the forces in snapshot by adding in the force bias from each
		// CV to each atom based on the gradient of the CV.
		for (size_t j = 0; j < forces.size(); ++j)
			{
			for(size_t k = 0; k < forces[j].size(); ++k)
				{
				if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[8] > 0 && _biases[j][k] != 0)
					_stringout << _biases[j][k] << " ";
				forces[j][k] += _biases[j][k];					
				}
			}
		
		// Feedback on locations.
		if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[6] > 0)
			{
			_stringout << " Locations are: ";
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
							_stringout << loc[j][k] << " ";
							}
						}
					}
				}
			}
			
		
	// Each output gets its own line
	if(snapshot->GetIteration() % _printdetails[0] == 0)
		_stringout << std::endl;	
	}

	// Post-simulation hook.
	void ABF::PostSimulation(Snapshot*, const CVList&)
	{

	_stringout << "End of simulation. The coordinates for Generalized Force Vector Field and N, in order of CVs in columns are:";

		int modulo = 1;
		int index;
		for(size_t i = 0; i < _N.size(); ++i)
			{
			index = i;
			_stringout << std::endl;
			for(size_t j=0 ; j < _histdetails.size(); ++j)
				{
				modulo = 1;
				for(size_t k=j+1 ; k <_histdetails.size(); ++k)
					{
					modulo = modulo * _histdetails[k][2];
					}
				_stringout << (floor(index/modulo)+0.5)*((_histdetails[j][1]-_histdetails[j][0])/_histdetails[j][2]) + _histdetails[j][0] << " ";
				index = index % modulo;
				}
			}
		_stringout << std::endl;
		_stringout << "Printing out 1) Dimensionality of F 2) Total Length of N 3) Total Length of F, 4) CV details 5) Initial histogrammed F estimate." << std::endl;
		_stringout << _dim << std::endl;
		_stringout << _Nworld.size() <<std::endl;
		_stringout << _Fworld.size() <<std::endl;
		_stringout << "CV bounds and nr of bins are: ";
		for(size_t i = 0; i < _histdetails.size(); ++i)
			{
			_stringout << "[" << _histdetails[i][0] << "," << _histdetails[i][1] << "," << _histdetails[i][2] << "] ";
			}
		_stringout << std::endl;
		for(size_t i = 0; i < _Nworld.size(); ++i)
				{
				for(size_t j = 0; j < _dim; ++j)
					{
					_stringout << _Fworld[_dim*i+j] << " ";
					}
				_stringout << _Nworld[i] << std::endl;
				}
			_stringout << std::endl;	
		_stringout.close();
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
				//_stringout << " Gradient[" << i << "]: ";
				for(size_t j = 0; j < _biases.size(); ++j)
					{
					for(size_t k = 0; k < _biases[j].size(); ++k)
						{
						if (grad[j][k] != 0)
							{
							//_stringout << grad[j][k] << " ";
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
						_stringout << " Minimum Harmonic Restraint Active on CV: " << i << "; Restraint force is: ";
					auto& grad = cvs[i]->GetGradient();
					for(size_t j = 0; j < _biases.size(); ++j)
					{
						for(size_t k = 0; k < _biases[j].size(); ++k)
						{
							if (grad[j][k] != 0)
							{
								_biases[j][k] += grad[j][k]*(_restraint[i][0] - cvs[i]->GetValue())*_restraint[i][2];
								if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[7] > 0)
									_stringout << _biases[j][k] << " ";
							}
						}
					}

				}

				else if(cvs[i]->GetValue() > _restraint[i][1] && _restraint[i][2] > 0)
				{
					if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[7] > 0)
						_stringout << " Maximum Harmonic Restraint Active on CV: " << i << "; Restraint force is: ";

					auto& grad = cvs[i]->GetGradient();
					for(size_t j = 0; j < _biases.size(); ++j)
					{
						for(size_t k = 0; k < _biases[j].size(); ++k)
						{
							if (grad[j][k] != 0)
							{
								_biases[j][k] += grad[j][k]*(_restraint[i][1] - cvs[i]->GetValue())*_restraint[i][2];
								if(snapshot->GetIteration() % _printdetails[0] == 0 && _printdetails[7] > 0)
									_stringout << _biases[j][k] << " ";
							}
						}
					}	
				}
			}
		}
	}
}



































