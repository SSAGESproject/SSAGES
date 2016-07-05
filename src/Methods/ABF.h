#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <vector>
#include <string>


namespace SSAGES
{

	//! Adaptive Biasing Force Algorithm
	/*!
	 * \ingroup Methods
	 *
	 *
	 * Implementation of the Adaptive Biasing Force algorithm based on:
	 * Darve, E.; Rodriguez-Gomez, D.; Pohorille, A. Adaptive biasing force method for scalar and vector free energy calculations, J. Chem. Phys. 2008, 128, 144120
	*/

	class ABF : public Method
	{
	private:	
		//! To store running total. 
		/*!
		 *A 1D vector, but will hold N-dimensional data, where N is number of CVs +1. This will be size (CVbinNr1*CVbinNr2*..)*3.
		 */
		std::vector<double> _F;

		//! Will hold the global total, synced to every time step. 
		/*!
		 *A 1D vector, but will hold N-dimensional data, where N is number of CVs +1. This will be size (CVbinNr1*CVbinNr2*..)*3.
		 */
		std::vector<double> _Fworld;

		//! To store number of hits at a given CV bin.
		/*!
		 *A 1D vector, but will hold N-dimensional data, where N is number of CVs. This will be size (CVbinNr1*CVbinNr2*..).
		 */
		std::vector<int> _N;

		//! To store number of hits at a given CV bin.
		/*!
		 *A 1D vector, but will hold N-dimensional data, where N is number of CVs. This will be size (CVbinNr1*CVbinNr2*..).
		 */
		std::vector<int> _Nworld;

		//! Information for a harmonic restraint to keep CV in the region of interest. 
		/*!
		 * This is a 2 Dimensional object set up as the following:
		 * _restraint is a vector of three vectors, each of those three vectors are (Number of CVs) long.
		 * First of these vectors hold the lower bound for the CV restraints in order.
		 * Second vector holds the upper bound for the CV restraints in order.
		 * Third vector holds the spring constants for the CV restraints in order.
		 */
		std::vector<std::vector<double>> _restraint;		

		//! The minimum number of hits required before full biasing, bias is _F[i]/max(_N[i],_min).
		int _min;

		//! To hold last iterations wdotp value for derivative
		std::vector<double> _wdotpold;

		//! To hold last iterations _F value for removing bias
		std::vector<double> _Fold;

		//! Function to return bin coordinate to address _F and _N, given a vector [CV1,CV2..] values
		int histCoords(const CVList& cvs);

		//! Thermodynamic beta.
		double _beta;

		//! To read in estimate of F.
		std::string _readF; 

		//! Biases.	
		std::vector<std::vector<double>> _biases;

		//! Number of CVs in system
		int _dim;

		//! Output stream for string data.
		std::ofstream _stringout;

		//! The node this belongs to
		unsigned int _mpiid;

		//! Histogram details. 
		/*!
		 * This is a 2 Dimensional object set up as the following:
		 * Histdetails is a vector of three vectors, each of those three vectors are (Number of CVs) long.
		 * First of these vectors hold the lower bound for the CVs in order.
		 * Second vector holds the upper bound for the CVs in order.
		 * Third vector holds number of histogram bins for the CVs in order.
		*/ 
		std::vector<std::vector<double>> _histdetails;

		//! Vector to hold print out information
		/*! 
		 * [Print every how many timesteps?, Print CVs?, Print Orthogonalization Correction?, Print Normalization Factor?, Print Gradient?, Print Genforce?, Print Coords?, Print Restraint info?, Print Biases?] 
		 */
		std::vector<int> _printdetails;

		//! Integer to hold F estimate backup information
		/*!
		 * Print every how many timesteps? ; -1: Do not backup during simulation
		 */
		int _FBackupInterv;

		//! Unit Conversion Constant from W dot P to Force
		/*!
		 * It is crucial that unit conversion is entered correctly.
                 * Unit conversion from d(momentum)/d(time) to force for the simulation. 
		 * For LAMMPS using units real, this is 2390.06 (gram.angstrom/mole.femtosecond^2 -> kcal/mole.angstrom)
		 */
		double _unitconv;
	
		//! Enable or Disable Gram-Schmidt Orthogonalization
		int _Orthogonalization;

		//! Computes the bias force.
		void CalcBiasForce(const CVList& cvs, const std::vector<double>& genforce, const Snapshot* snapshot);
		
		//! Timestep of integration
		double _timestep;

	public: 
		//! Constructor 
		/*!
		 * Constructs an instance of Adaptive Biasing Force method.
		 * _histdetails holds three arrays that define the minimum, the maximum, and the number of bins for each CV.
		 * _restraint holds the minimum, the maximum and the spring constant for CV restraints. These restraints should be outside of the range defined in _histdetails by at least one bin size on each side.
		 * _timestep is the time between two steps in the simulations, in simulation units.
		 * _min is the minimum number of hits in a histogram bin required before full biasing is applied.
		 * _readF is used to provide a prior histogram to restart a simulation, or to provide a guess.
		 * _printdetails is used to set up what information to print and the frequency of printing out.
		 * _FBackupInterv is used to set up how often the adaptive force histogram is saved.
		 * _unitconv is used to provide the unit conversion from d(momentum)/d(time) to force.
		 * _Orthogonalization is a flag to turn on or off Gram-Schmidt Orthogonalization.
		 */ 

		ABF(boost::mpi::communicator& world,
			 boost::mpi::communicator& comm,
			 const std::vector<std::vector<double>>& histdetails, std::vector<std::vector<double>> restraint, double timestep, double min, std::string readF, std::vector<int> printdetails, int FBackupInterv, double unitconv, int Orthogonalization, unsigned int frequency) : 
		Method(frequency, world, comm), _biases(0), _wdotpold(0), _mpiid(0), _Fold(0), _beta(0), _dim(0), _F(0), _N(0), _Fworld(0), _Nworld(0), _histdetails(histdetails), _restraint(restraint), _timestep(timestep), _min(min), _readF(readF), _printdetails(printdetails), _FBackupInterv(FBackupInterv), _unitconv(unitconv), _Orthogonalization(Orthogonalization)
		{
		}

		//! Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		void Serialize(Json::Value& json) const override
		{
		
		}
		//! Destructor
		~ABF() {}
	};
}
			
