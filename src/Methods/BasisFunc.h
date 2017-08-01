/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Joshua Moller <jmoller@uchicago.edu>
 *           2017 Julian Helfferich <julian.helfferich@gmail.com>
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
#include "Utility/Basis.h"
#include "Grids/Histogram.h"
#include <fstream>
#include <vector>

namespace SSAGES
{		
    //! Basis Function Sampling Algorithm
    /*!
     * \ingroup Methods
     *
     * Implementation of the Basis Function Sampling Method based on
     * \cite WHITMER2014190602.
     */
	class BFS : public Method
	{
	private:	
        
        //! Histogram of visited states.
        /*!
         * Histogram is stored locally.
         */
        Grid<uint> *h_;

        //! Stored bias potential
        /*!
         * The sum of basis functions that adds up to the bias potential of the surface
         */
        Grid<double> *b_;

        //! Stored gradients
        /*!
         * The gradients corresponding to each point on the grid
         */
        Grid<std::vector<double>> *f_;

        //! The basis evaluator class
        /*!
         * A class which holds the basis functions and updates the coefficients
         * at each cycle
         */
        BasisEvaluator evaluator_;

        //! The biased histogram of states.
        /*!
         * The biased histogram of states has the form hist_*exp(phi*beta),
         * where phi is the bias potential and beta is the inverse of the
         * temperature. It is defined globally.
         */
        std::vector<double> unbias_;

        //! The coefficient array for restart runs
        std::vector<double> coeff_arr_;

        //! Spring constants for restrained system.
        /*!
         * The system uses this to determine if the system is to be restrained
         * on the defined interval. The user inputs the spring constants if the
         * system is not periodic.
         */
        std::vector<double> restraint_;

        //! Upper position of the spring restraint.
        std::vector<double> boundUp_;

        //! Lower position of the spring restraint.
        std::vector<double> boundLow_;
        
        //! Frequency of coefficient updates
		unsigned int cyclefreq_;
        
        //! The node that the current system belongs to, primarily for printing and debugging.
        unsigned int mpiid_;

        //! Weighting for potentially faster sampling.
        /*!
         * Weighting can be used to potentially sample faster, however it can
         * cause the simulation to explode. By default this value will be set
         * to 1.0
         */
        double weight_;

        //! Self-defined temperature.
        /*!
         * In the case of the MD engine using a poorly defined temperature, the
         * option to throw it into the method is available. If not provided it
         * takes the value from the engine.
         */
        double temperature_;

        //! The tolerance criteria for the system .
        double tol_;

        //! A variable to check to see if the simulation is in bounds or not.
        bool bounds_;

        //! A check to see if you want the system to end when it reaches the convergence criteria.
        bool converge_exit_;

        //! The functions which calculates the updated bias and coefficients and then stores them
        void UpdateBias(const CVList& cvs, const double beta);

        //! A function which checks to see if the CVs are still in bounds
        void InBounds(const CVList& cvs);

		//! Prints the current bias to a defined file from the JSON.
        /*!
         * \param cvs List of collective variables.
         * \param beta Scale parameter.
         */
		void PrintBias(const CVList& cvs, const double beta);

		//! Output stream for basis projection data.
		std::ofstream basisout_;

        //! The option to name both the basis and coefficient files will be given
        //! Basis filename 
        std::string bnme_;

        //! Iteration counter. 
        uint iteration_;

	public:
        //! Constructor
		/*!
         * \param world MPI global communicator.
         * \param comm MPI local communicator.
         * \param polyord Order of Legendre polynomials.
         * \param restraint Restraint spring constants.
         * \param boundUp Upper bounds of restraint springs.
         * \param boundLow Lower bounds of restraint springs.
         * \param cyclefreq Cycle frequency.
         * \param frequency Frequency with which this Method is applied.
         * \param bnme Basis file name.
         * \param cnme Coefficient file name.
         * \param temperature Automatically set temperature.
         * \param tol Threshold for tolerance criterion.
         * \param weight Weight for improved sampling.
         * \param converge If \c True quit on convergence.
         *
         * Constructs an instance of the Basis function sampling method. The
         * coefficients describes the basis projection of the system. This is
         * updated once every cyclefreq_. For now, only the Legendre polynomial
         * is implemented. Others will be added later.
         */
		BFS(const MPI_Comm& world,
			  const MPI_Comm& comm,
              Grid<uint> *h,
              Grid<std::vector<double>> *f,
              Grid<double> *b,
              const std::vector<BasisFunction>& functions,
              const std::vector<double>& restraint,
              const std::vector<double>& boundUp,
              const std::vector<double>& boundLow,
              unsigned int cyclefreq,
			  unsigned int frequency,
              const std::string bnme,
              const double temperature,
              const double tol,
              const double weight,
              bool converge) :
		Method(frequency, world, comm), 
        h_(h),  b_(b), f_(f), unbias_(), coeff_arr_(), evaluator_(functions),
        restraint_(restraint), boundUp_(boundUp), boundLow_(boundLow),
        cyclefreq_(cyclefreq), mpiid_(0), weight_(weight),
        temperature_(temperature), tol_(tol),
        converge_exit_(converge), bnme_(bnme), iteration_(0)
		{
		}

		//! Pre-simulation hook.
        /*!
         * \param snapshot Simulation snapshot.
         * \param cvmanager Collective variable manager.
         */
		void PreSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Post-integration hook.
        /*!
         * \param snapshot Simulation snapshot.
         * \param cvmanager Collective variable manager.
         */
		void PostIntegration(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Post-simulation hook.
        /*!
         * \param snapshot Simulation snapshot.
         * \param cvmanager Collective variable manager.
         */
		void PostSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

        //! Set the current iteration
        /*!
         * \param iter New value for the current iteration.
         *
         * This function is used to set the current iteration, for example when
         * continuing from a restart.
         */
        void SetIteration(const int iter)
        {
            iteration_ = iter;
        }

        //! Set the values for the basis.
        /*!
         * \param coeff List of coefficients.
         * \param unbias List of unbiased values.
         *
         * This function is used to set starting values at the beginning of
         * a run. For example when continuing from a restart value.
         */
        void SetBasis(const std::vector<double>&coeff, std::vector<double>&unbias)
        {
            coeff_arr_ = coeff;
            unbias_ = unbias;
        }

		//! \copydoc Buildable::Build()
		static BFS* Build(const Json::Value& json, 
		                        const MPI_Comm& world,
		                        const MPI_Comm& comm,
					            const std::string& path);

        //! Destructor.
        ~BFS()
        {
            delete h_;
            delete f_;
            delete b_;
        }
	};
}
			
