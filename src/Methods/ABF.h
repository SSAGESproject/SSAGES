#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <vector>
#include <string>

#include <Eigen/Dense>

namespace SSAGES
{

	//! Adaptive Biasing Force Algorithm
	/*!
	 * \ingroup Methods
	 *
	 *
	 * Implementation of the Adaptive Biasing Force algorithm based on
	 * \cite DARVE2008144120
	 */

	class ABF : public Method
	{
	private:	
		//! To store running total. 
		/*!
		 * A 1D vector, but will hold N-dimensional data, where N is number of
		 * CVs +1. This will be size (CVbinNr1*CVbinNr2*..)*3.
		 */
		Eigen::VectorXd _F;

		//! Will hold the global total, synced to every time step. 
		/*!
		 * A 1D vector, but will hold N-dimensional data, where N is number of
		 * CVs +1. This will be size (CVbinNr1*CVbinNr2*..)*3.
		 */
		Eigen::VectorXd _Fworld;

		//! To store number of hits at a given CV bin.
		/*!
		 * A 1D vector, but will hold N-dimensional data, where N is number of
		 * CVs. This will be size (CVbinNr1*CVbinNr2*..).
		 */
		std::vector<int> _N;

		//! To store number of hits at a given CV bin.
		/*!
		 * A 1D vector, but will hold N-dimensional data, where N is number of
		 * CVs. This will be size (CVbinNr1*CVbinNr2*..).
		 */
		std::vector<int> _Nworld;

		//! Information for a harmonic restraint to keep CV in the region of interest. 
		/*!
		 * This is a 2 Dimensional object set up as the following:
		 * _restraint is a vector of (nr of CV) vectors, ordered in CV order.
		 * Each of those vectors are 3 long.
		 * First entry of each vector holds the lower bound for that CV restraint.
		 * Second entry of each vector holds the upper bound for that CV restraint.
		 * Third entry of each vector holds the spring constant for that CV restraint.
		 */
		std::vector<std::vector<double>> _restraint;		

		//! The minimum number of hits required before full biasing, bias is
		//!_F[i]/max(_N[i],_min).
		int _min;

		//! To hold last two iterations wdotp value for derivative
		Eigen::VectorXd _wdotp1, _wdotp2;

		//! To hold last iterations _F value for removing bias
		Eigen::VectorXd _Fold;

		//! Get coordinates of histogram bin corresponding to given list of CVs.
		/*!
		 * \param cvs List of CVs.
		 *
		 * \return Index of histogram bin.
		 *
		 * Function to return bin coordinate to address _F and _N, given a vector
		 * [CV1,CV2..] values.
		 */
		int histCoords(const CVList& cvs);

		//! Thermodynamic beta.
		double _beta;

		//! Biases.	
		std::vector<Vector3> _biases;

		//! Number of CVs in system
		unsigned int _dim;

		//! Output stream for walker-specific data.
		std::ofstream _walkerout;

		//! Output stream for world data.
		std::ofstream _worldout;

		//! File name for world data
		std::string _filename;

		//! The node this belongs to
		unsigned int _mpiid;

		//! Histogram details. 
		/*!
		 * This is a 2 Dimensional object set up as the following:
		 * Histdetails is a vector of (nr of CV) vectors, ordered in CV order.
		 * Each of those vectors are 3 long.
		 * First entry of each vector holds the lower bound for that CV.
		 * Second entry of each vector holds the upper bound for that CV.
		 * Third entry of each vector holds the nr of bins for that CV dimension.
		*/ 
		std::vector<std::vector<double>> _histdetails;

		//! Vector to hold print out information
		/*! 
		 * [Print every how many timesteps?, Print CVs?, Print
		 * Orthogonalization Correction?, Print Normalization Factor?, Print
		 * Gradient?, Print Genforce?, Print Coords?, Print Restraint info?,
		 * Print Biases?]
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
		 * For LAMMPS using units real, this is
		 * 2390.06 (gram.angstrom/mole.femtosecond^2 -> kcal/mole.angstrom)
		 */
		double _unitconv;
	
		//! Enable or Disable Gram-Schmidt Orthogonalization
		int _Orthogonalization;

		//! Computes the bias force.
		void CalcBiasForce(const Snapshot* snapshot, const CVList& cvs, int coord);
		
		//! Writes out data to file.
		void WriteData();

		//! Timestep of integration
		double _timestep;

	public: 
		//! Constructor 
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param histdetails Minimum, maximum and number of bins for each CV.
		 * \param restraint Minimum, maximum and spring constant for CV restraints.
		 * \param timestep Simulation time step.
		 * \param min Minimum number of hist in a histogram bin before biasing is applied.
		 * \param filename Name for output file.
		 * \param printdetails Set up what information to print and frequency of printing.
		 * \param FBackupInterv Set how often the adaptive force histogram is saved.
		 * \param unitconv Unit conversion from d(momentum)/d(time) to force.
		 * \param Orthogonalization Flag to turn on or off Gram-Schmidt Orthogonalization.
		 * \param frequency Frequency with which this method is invoked.
		 *
		 * Constructs an instance of Adaptive Biasing Force method.
		 * \note The restraints should be outside of the range defined in
		 *       _histdetails by at least one bin size on each side.
		 */ 
		ABF(boost::mpi::communicator& world,
			boost::mpi::communicator& comm,
			const std::vector<std::vector<double>>& histdetails,
			std::vector<std::vector<double>> restraint,
			double timestep,
			double min,
			std::string filename,
			std::vector<int> printdetails,
			int FBackupInterv,
			double unitconv,
			int Orthogonalization,
			unsigned int frequency) :
		Method(frequency, world, comm), _F(), _Fworld(), _N(0), _Nworld(0),
		_restraint(restraint), _min(min), _wdotp1(), _wdotp2(), _Fold(), _beta(0),
		_filename(filename), _biases(), _dim(0), _mpiid(0), _histdetails(histdetails),
		_printdetails(printdetails), _FBackupInterv(FBackupInterv),
		_unitconv(unitconv), _Orthogonalization(Orthogonalization),
		_timestep(timestep)
		{
		}

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Set biasing histogram
		/*!
		 * \param F Vector containing values for the running total.
		 * \param N Vector containing number of hits for bin intervals.
		 */
		void SetHistogram(const Eigen::VectorXd& F, const std::vector<int>& N)
		{
			_F = F;
			_N = N;
		}		
		
		//! Set iteration of the method
		/*!
		 * \param iter New value for the iteration counter.
		 */
		void SetIteration(const int iter)
		{
			_iteration = iter;
		}			

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override
		{
		
			json["type"] = "ABF";

			for(size_t i = 0; i < _histdetails.size(); ++i)
			{
				json["CV_lower_bounds"].append(_histdetails[i][0]);				
				json["CV_upper_bounds"].append(_histdetails[i][1]);
				json["CV_bins"].append(_histdetails[i][2]);
			}

			for(size_t i = 0; i < _restraint.size(); ++i)
			{
				json["CV_restraint_minimums"].append(_restraint[i][0]);
				json["CV_restraint_maximums"].append(_restraint[i][1]);
				json["CV_restraint_spring_constants"].append(_restraint[i][2]);
			}

			json["timestep"] = _timestep;

			json["minimum_count"] = _min;

			for(size_t i = 0; i < _printdetails.size(); ++i)
				json["print_details"].append(_printdetails[i]);
			 
			json["backup_frequency"] = _FBackupInterv;
			
			json["unit_conversion"] = _unitconv;

			json["orthogonalization"] = _Orthogonalization;			
		
			for(size_t i = 0; i < _F.size(); ++i)
				json["F"].append(_F[i]);

			for(size_t i = 0; i < _N.size(); ++i)
				json["N"].append(_N[i]);

			json["iteration"] = _iteration;
			
			json["filename"] = _filename;		
		
		}

		//! Destructor
		~ABF() {}
	};
}
			
