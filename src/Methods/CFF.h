/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2018 Emre Sevgen <sesevgen@uchicago.edu> and Hythem Sidky <hsidky@nd.edu>
 *           2020 Elizabeth M.Y. Lee <emlee@uchicago.edu> and Boyuan Yu <boyuanyu@uchicago.edu>
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
#include "Grids/Grid.h"
#include "nnet/nnet.h"

namespace SSAGES
{
	//! Combined Force Frequency (CFF) Algorithm.
	/*!
	 * \ingroup Methods
	 *
	 * Implementation of the CFF algorithm based on
	 * Sevgen et al. J. Chem. Theory Comput. 2020, 16, 3, 1448-1455.
	 */
	class CFF : public Method
	{
	private:
		//! Neural network topology. 
		Eigen::VectorXi topol_;

		//! Number of iterations per sweep. 
		unsigned int sweep_, nsweep_;

		//! Number of iterations after which we turn on full weight. 
		unsigned int citers_; 

		//! Neural network trained using visit frequency
		nnet::neural_net net_; 
        
        	//! Neural network trained on both visit frequency and force
		nnet::neural_net net2_;

		//! Timestep.
		double timestep_;

		//! Previous and current histogram weight. 
		double pweight_, weight_;

		//! System temperature and energy units. 
		double temp_, kbt_;

		//! Generalized force grid that stores total of the local walker.
		std::vector<Grid<double>*> F_;

       	//! Generalized force grid that stors the global total.
		std::vector<Grid<double>*> Fworld_;
	
		//! To hold the last iterations F_ value for removing bias.
		Eigen::VectorXd Fold_;
   
       	//! To hold last iteration wdotp value for numerical derivative.
        Eigen::VectorXd wdotp1_; 

        //! To hold second to last iteration wdotp value for numerical derivative.
        Eigen::VectorXd wdotp2_;

		//! Force grid. 
		Grid<Eigen::VectorXd>* fgrid_;

		//! Histogram grid.
		Grid<unsigned int>* hgrid_;

		//! Histogram grid for that sotres the number of global hits
		Eigen::ArrayXi Nworld_;

		//! Unbiased histogram grid. 
		Grid<double>* ugrid_;

		//! Hold parameters to adjust ratio of 1st and 2nd neural networks (freq vs force-based).
		double ratio_;

		//! Eigen matrices of grids. 
		Eigen::MatrixXd hist_, bias_;

		//! Bounds 
		std::vector<double> lowerb_, upperb_;

		//! Bound restraints. 
		std::vector<double> lowerk_, upperk_;

		//! Output filename. 
		std::string outfile_;

		//! Overwrite outputs. 
		bool overwrite_;

		//! Trains the neural network.
		void TrainNetwork();

		//! Writes out the bias and CFF restart data to file. 
		void WriteBias();

		//! Number of CVs in system.
		int dim_;

		//! Unit conversion from mass*velocity/time to force.
		double unitconv_;

        //! The minimum number of hits required before full biasing, bias is
        //! F_[i]/max(N_[i],min_).
        int min_;

        //! To hold booleans for training neural nerwork only in specific region for net2_.
        Eigen::MatrixXd force_to_val_ratio_;

        //! Restart from CFF.
        bool restart_from_cff_;

        //! Restart from ABF.
        bool restart_from_abf_;

	public: 
		//! Constructor 
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param topol Topology of network.
		 * \param fgrid Grid containing biasing forces.
		 * \param hgrid Grid containing histogram.
		 * \param ugrid Grid containing unbiased histogram.
		 * \param lowerb Lower bounds for CVs.
		 * \param upperb Upper bounds for CVs.
		 * \param lowerk Lower bound restraints for CVs.
		 * \param upperk Upper bound restraints for CVs.
		 * \param temperature Temperature of the simulation.
		 * \param weight Relative weight of the statistics in sweep.
		 * \param nsweep Number of iterations in the sweep.
		 * \param min Minimum number of hits required for the full force biasing
         	 * \param restart_from_cff Restart flag starting from previous CFF runs
         	 * \param restart_from_abf Restart flag starting from previous ABF runs
		 *
		 * Constructs an instance of Combined Force Frequency method.
		 */ 
		CFF(const MPI_Comm& world, 
		    const MPI_Comm& comm, 
		    const Eigen::VectorXi& topol,
		    Grid<Eigen::VectorXd>* fgrid,
		    Grid<unsigned int>* hgrid,
		    Grid<double>* ugrid,
		    std::vector<Grid<double>*> F, std::vector<Grid<double>*> Fworld,
		    const std::vector<double>& lowerb,
		    const std::vector<double>& upperb,
		    const std::vector<double>& lowerk,
		    const std::vector<double>& upperk,
		    double temperature,
		    double unitconv,
		    double timestep,
		    double weight,
		    unsigned int nsweep,
		    int min,
		    bool restart_from_cff,
		    bool restart_from_abf
		    );

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		void PreSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		void PostIntegration(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		void PostSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Set previous history weight.
		void SetPrevWeight(double h)
		{
			pweight_ = h;
		}

		//! Set name of output file. 
		void SetOutput(const std::string& outfile)
		{
			outfile_ = outfile;
		}

		//! Set overwrite flag on output file.
		void SetOutputOverwrite(bool overwrite)
		{
			overwrite_ = overwrite;
		}

		//! Set number of iterations after which we turn on full weight. 
		void SetConvergeIters(unsigned int citers)
		{
			citers_ = citers;
		}

		//! Set maximum number of training iterations per sweep.
		void SetMaxIters(unsigned int iters)
		{
			auto params = net_.get_train_params();
			params.max_iter = iters;
			net_.set_train_params(params);

			auto params2 = net2_.get_train_params();
			params2.max_iter = iters;
			net2_.set_train_params(params2);
		}

		//! Set minimum loss function value (should be zero for production).
		void SetMinLoss(double loss)
		{
			auto params = net_.get_train_params();
			params.min_loss = loss; 
			net_.set_train_params(params);
		
			auto params2 = net2_.get_train_params();
			params2.min_loss = loss; 
			net2_.set_train_params(params2);			
		}
		
		//! Set training ratio for gradient vs value.
		void SetRatio(double trainratio)
		{
			auto params = net_.get_train_params();
			params.ratio = trainratio; 
			net_.set_train_params(params);

			auto params2 = net2_.get_train_params();
			params2.ratio = trainratio; 
			net2_.set_train_params(params2);		
		}

		//! \copydoc Method::BuildMethod().
		static CFF* Build(
			const Json::Value& json,
			const MPI_Comm& world,
			const MPI_Comm& comm,
			const std::string& path);
		
        	//! Destructor.
		~CFF()
		{
			delete fgrid_; 
			delete hgrid_;
		}
	};
}
