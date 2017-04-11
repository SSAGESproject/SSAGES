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
#include "CVs/CollectiveVariable.h"
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
	 * Implementation of the Adaptive Biasing Force algorithm based on
	 * \cite DARVE2008144120
	 */
	class ABF : public Method, public Buildable<ABF>
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
		Eigen::VectorXd Fworld_;

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
		std::vector<int> Nworld_;

		//! Information for a harmonic restraint to keep CV in the region of interest. 
		/*!
		 * This is a 2 Dimensional object set up as the following:
		 * restraint_ is a vector of (nr of CV) vectors, ordered in CV order.
		 * Each of those vectors are 3 long.
		 * First entry of each vector holds the lower bound for that CV restraint.
		 * Second entry of each vector holds the upper bound for that CV restraint.
		 * Third entry of each vector holds the spring constant for that CV restraint.
		 */
		std::vector<std::vector<double>> restraint_;

		//! For each CV, holds whether that CV has periodic boundaries or not.
		std::vector<bool> isperiodic_;

		//! Holds periodic boundaries of CVs.
		/*!
		 * This is a 2 Dimensional object set up as the following:
		 * periodicboundaries_ is a vector of (nr of CV) vectors, ordered in CV order.
		 * Each of those vectors are 2 long.
		 * First entry of each vector holds the lower bound for that CV periodic boundary.
		 * Second entry of each vector holds the upper bound for that CV periodic boundary.
		 */
		std::vector<std::vector<double>> periodicboundaries_;		

		//! The minimum number of hits required before full biasing, bias is
		//!_F[i]/max(_N[i],min_).
		int min_;

		//! To hold last two iterations wdotp value for derivative
		Eigen::VectorXd wdotp1_, wdotp2_;

		//! To hold last iterations _F value for removing bias
		Eigen::VectorXd Fold_;

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

		//! Mass weighing of bias enabled/disabled
		bool massweigh_;

		//! Biases.	
		std::vector<Vector3> biases_;

		//! Number of CVs in system
		unsigned int dim_;

		//! Output stream for world data.
		std::ofstream worldout_;

		//! File name for world data
		std::string filename_;

		//! Histogram details. 
		/*!
		 * This is a 2 Dimensional object set up as the following:
		 * Histdetails is a vector of (nr of CV) vectors, ordered in CV order.
		 * Each of those vectors are 3 long.
		 * First entry of each vector holds the lower bound for that CV.
		 * Second entry of each vector holds the upper bound for that CV.
		 * Third entry of each vector holds the nr of bins for that CV dimension.
		*/ 
		std::vector<std::vector<double>> histdetails_;

		//! Integer to hold F estimate backup information
		/*!
		 * Print every how many timesteps? ; -1: Do not backup during simulation
		 */
		int FBackupInterv_;

		//! Unit Conversion Constant from W dot P to Force
		/*!
		 * It is crucial that unit conversion is entered correctly.
		 * Unit conversion from d(momentum)/d(time) to force for the simulation.
		 * For LAMMPS using units real, this is
		 * 2390.06 (gram.angstrom/mole.femtosecond^2 -> kcal/mole.angstrom)
		 */
		double unitconv_;

		//! Computes the bias force.
		void CalcBiasForce(const Snapshot* snapshot, const CVList& cvs, int coord);
		
		//! Writes out data to file.
		void WriteData();

		//! Timestep of integration
		double timestep_;

		//! Iteration counter. 
		uint iteration_;

		//! Mass vector. Empty unless required.
		Eigen::VectorXd mass_;

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
		 *       histdetails_ by at least one bin size on each side.
		 */ 
		ABF(const MPI_Comm& world,
			const MPI_Comm& comm,
			std::vector<std::vector<double>> restraint,
			std::vector<bool> isperiodic,
			std::vector<std::vector<double>> periodicboundaries,
			double min,
			bool massweigh,
			std::string filename,
			const std::vector<std::vector<double>>& histdetails,
			int FBackupInterv,
			double unitconv,
			double timestep,
			unsigned int frequency) :
		Method(frequency, world, comm), _F(), Fworld_(), _N(0), Nworld_(0),
		restraint_(restraint), isperiodic_(isperiodic), periodicboundaries_(periodicboundaries),
		min_(min), wdotp1_(), wdotp2_(), Fold_(), massweigh_(massweigh),
		biases_(), dim_(0), filename_(filename), histdetails_(histdetails), 
		FBackupInterv_(FBackupInterv), unitconv_(unitconv), timestep_(timestep),
		iteration_(0)
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
			iteration_ = iter;
		}			
		
		//! \copydoc Buildable::Build()
		static ABF* Build(const Json::Value& json, 
		                  const MPI_Comm& world,
		                  const MPI_Comm& comm,
					      const std::string& path);

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override;

		//! Destructor
		~ABF() {}
	};
}
			
