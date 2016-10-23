/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Emre Sevgen <sesevgen@uchicago.edu>
 *                Joshua Moller <jmoller@uchicago.edu>
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

		//! To hold the last iteration wdotp value for derivative.
		Eigen::VectorXd _wdotp1

		//! To hold the one-but-last iteration wdotp value for derivative.
		Eigen::VextorXd _wdotp2;

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

		//! Computes the bias force.
		void CalcBiasForce(const Snapshot* snapshot, const CVList& cvs, int coord);
		
		//! Writes out data to file.
		void WriteData();

		//! Timestep of integration
		double _timestep;

		//! Number of cvs
		int _ncv = 0;

		//! Mass vector. Empty unless required.
		Eigen::VectorXd _mass;

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
		 * \param FBackupInterv Set how often the adaptive force histogram is saved.
		 * \param unitconv Unit conversion from d(momentum)/d(time) to force.
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
			int FBackupInterv,
			double unitconv,
			unsigned int frequency) :
		Method(frequency, world, comm), _F(), _Fworld(), _N(0), _Nworld(0),
		_restraint(restraint), _min(min), _wdotp1(), _wdotp2(), _Fold(), _beta(0),
		_biases(), _dim(0), _filename(filename), _mpiid(0), _histdetails(histdetails), 
		_FBackupInterv(FBackupInterv), _unitconv(unitconv), _timestep(timestep)
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
			json["backup_frequency"] = _FBackupInterv;			
			json["unit_conversion"] = _unitconv;
			json["iteration"] = _iteration;
			json["filename"] = _filename;		
				
			for(int i = 0; i < _F.size(); ++i)
				json["F"].append(_F[i]);

			for(size_t i = 0; i < _N.size(); ++i)
				json["N"].append(_N[i]);
	
		
		}

		//! Destructor
		~ABF() {}
	};
}
			
