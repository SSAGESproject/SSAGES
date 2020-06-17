/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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
	//! Artificial Neural Network Method
	/*!
	 * \ingroup Methods
	 *
	 * Implementation of the Artificial Neural Network Method based on
	 * \cite SIDKY2018104111
	 */
	class ANN : public Method
	{
	private:
		//! Neural network topology.
		Eigen::VectorXi topol_;

		//!@{
		//! Number of iterations per sweep.
		unsigned int sweep_, nsweep_;
		//!@}

		//! Number of iterations after which we turn on full weight.
		unsigned int citers_;

		//! Neural network.
		nnet::neural_net net_;

		//!@{
		//! Previous and current histogram weight.
		double pweight_, weight_;
		//!@}

		//!@{
		//! System temperature and energy units.
		double temp_, kbt_;
		//!@}

		//! Force grid.
		Grid<Eigen::VectorXd>* fgrid_;

		//! Histogram grid.
		Grid<unsigned int>* hgrid_;

		//! Unbiased histogram grid.
		Grid<double>* ugrid_;

		//!@{
		//! Eigen matrices of grids.
		Eigen::MatrixXd hist_, bias_, rbias_;
		//!@}

		//!@{
		//! Bounds
		std::vector<double> lowerb_, upperb_;
		//!@}

		//!@{
		//! Bound restraints.
		std::vector<double> lowerk_, upperk_;
		//!@}

		//! Output filename.
		std::string outfile_;

		//! Is the network preloaded?
		bool preloaded_;

		//! Overwrite outputs?
		bool overwrite_;

		//! Trains the neural network.
		void TrainNetwork();

		//! Writes out the bias to file.
		void WriteBias();

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
		 *
		 * Constructs an instance of Artificial Neural Network method.
		 */ 
		ANN(const MPI_Comm& world, 
		    const MPI_Comm& comm, 
		    const Eigen::VectorXi& topol,
		    Grid<Eigen::VectorXd>* fgrid,
		    Grid<unsigned int>* hgrid,
		    Grid<double>* ugrid,
		    const std::vector<double>& lowerb,
		    const std::vector<double>& upperb,
		    const std::vector<double>& lowerk,
		    const std::vector<double>& upperk,
		    double temperature,
		    double weight,
		    unsigned int nsweep
		);

		//! \copydoc Method::PreSimulation()
		void PreSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! \copydoc Method::PostIntegration()
		void PostIntegration(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! \copydoc Method::PostSimulation()
		void PostSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Set previous history weight.
		/*!
		 * \param h History weight
		 */
		void SetPrevWeight(double h)
		{
			pweight_ = h;
		}

		//! Set name of output file. 
		/*!
		 * \param outfile Output file
		 */
		void SetOutput(const std::string& outfile)
		{
			outfile_ = outfile;
		}

		//! Set overwrite flag on output file.
		/*!
		 * \param overwrite Boolean if output file should be overwritten
		 */
		void SetOutputOverwrite(bool overwrite)
		{
			overwrite_ = overwrite;
		}

		//! Set number of iterations after which we turn on full weight. 
		/*!
		 * \param citers Number of iterations before full weight
		 */
		void SetConvergeIters(unsigned int citers)
		{
			citers_ = citers;
		}

		//! Set maximum number of training iterations per sweep.
		/*!
		 * \param iters Maximum iterations per sweep
		 */
		void SetMaxIters(unsigned int iters)
		{
			auto params = net_.get_train_params();
			params.max_iter = iters;
			net_.set_train_params(params);
		}

		//! Set minimum loss function value (should be zero for production).
		/*!
		 * \param loss Minimum loss function value
		 */
		void SetMinLoss(double loss)
		{
			auto params = net_.get_train_params();
			params.min_loss = loss; 
			net_.set_train_params(params);
		}

		//! Load network state and bias from file.
		void ReadBias(const std::string&, const std::string&);

		//! \copydoc Method::BuildMethod()
		static ANN* Build(
			const Json::Value& json,
			const MPI_Comm& world,
			const MPI_Comm& comm,
			const std::string& path);

		~ANN()
		{
			delete fgrid_;
			delete hgrid_;
		}
	};
}
