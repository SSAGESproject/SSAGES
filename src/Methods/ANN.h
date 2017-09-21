#pragma once 

#include "Method.h"
#include "Grids/Grid.h"
#include "nnet/nnet.h"

namespace SSAGES
{
	class ANN : public Method
	{
	private:
		// Neural network topology. 
		Eigen::VectorXi topol_;

		//! Number of iterations per sweep. 
		uint sweep_, nsweep_;

		//! Neural network.
		nnet::neural_net net_;

		//! Previous and current histogram weight. 
		double pweight_, weight_;

		//! System temperature and energy units. 
		double temp_, kbt_;

		//! Force grid. 
		Grid<Eigen::VectorXd>* fgrid_;

		//! Histogram grid.
		Grid<uint>* hgrid_;

		//! Unbiased histogram grid. 
		Grid<double>* ugrid_;

		//! Eigen matrices of grids. 
		Eigen::MatrixXd hist_, bias_;
	   
		//! Bounds 
	   	std::vector<double> lowerb_, upperb_; 
	   
		//! Bound restraints. 
		std::vector<double> lowerk_, upperk_;

		//! Trains the neural network.
		void TrainNetwork();
		 
		//! Writes out the bias to file. 
		void WriteBias();

	public:
		ANN(const MPI_Comm& world, 
			const MPI_Comm& comm, 
			const Eigen::VectorXi& topol,
			Grid<Eigen::VectorXd>* fgrid,
			Grid<uint>* hgrid,
			Grid<double>* ugrid,
			const std::vector<double>& lowerb,
			const std::vector<double>& upperb,
			const std::vector<double>& lowerk,
			const std::vector<double>& upperk,
            double temperature,
            double weight,
            uint nsweep
		);

			//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		 void PreSimulation(Snapshot*, const class CVManager&) override;
		 
		//! Post-integration hook.
		/*!
		* \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		*/
		void PostIntegration(Snapshot*, const class CVManager&) override;

		//! Post-simulation hook.
		/*!
		* \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		*/
		void PostSimulation(Snapshot*, const class CVManager&) override;

		//! Set previous history weight.
        void SetPrevWeight(double h)
        {
            pweight_ = h;
		}

		//! \copydoc Buildable::Build()
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