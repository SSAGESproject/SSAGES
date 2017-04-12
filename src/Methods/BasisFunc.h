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
#include "Grids/Histogram.h"
#include <fstream>
#include <vector>

namespace SSAGES
{

    //! Map for histogram and coefficients.
    /*!
     * A clean mapping structure for both the histogram and the coefficients.
     * All vectors are written as 1D with a row major mapping. In order to make
     * iterating easier, the mapping of the 1D vectors are written here.
     */
    struct Map
    {
        //! The coefficient value
        double value;

        //! The mapping in an array of integers
        std::vector<int> map;

        //! Constructor
        /*!
         * \param map The mapping in an array of integers.
         * \param value The coefficient value.
         */
        Map(const std::vector<int>& map,
            double value) :
            value(value), map(map)
        {}
    };

    //! Look-up table for basis functions.
	/*!
     * The structure that holds the Look-up table for the basis function. To
     * prevent repeated calculations, both the derivatives and values of the
     * Legendre polynomials is stored here. More will be added in future
     * versions.
     */
	struct BasisLUT
	{
		//! The values of the basis sets
		std::vector<double> values;

		//! The values of the derivatives of the basis sets
		std::vector<double> derivs;

        //! Constructor.
        /*!
         * \param values The values of the basis sets.
         * \param derivs The values of the derivatives of the basis sets.
         */
		BasisLUT(const std::vector<double>& values,
			const std::vector<double>& derivs) :
			values(values), derivs(derivs)
		{}
	};
		
    //! Basis Function Sampling Algorithm
    /*!
     * \ingroup Methods
     *
     * Implementation of the Basis Function Sampling Method based on
     * \cite WHITMER2014190602.
     */
	class Basis : public Method, public Buildable<Basis>
	{
	private:	
        
        //! Histogram of visited states.
        /*!
         * Histogram is stored locally.
         */
        Histogram<int> *hist_;

        //! Globally located coefficient values.
		/*!
         * As coefficients are updated at the same time, the coefficients
         * should only be defined globally.
         */
		std::vector<Map> coeff_;

        //! The biased histogram of states.
        /*!
         * The biased histogram of states has the form hist_*exp(phi*beta),
         * where phi is the bias potential and beta is the inverse of the
         * temperature. It is defined globally.
         */
        std::vector<double> unbias_;

        //! The coefficient array for restart runs
        std::vector<double> coeff_arr_;

        //! The Basis set lookup table, also defined globally
		std::vector<BasisLUT> LUT_;

        //! Derivatives of the bias potential
        /*!
         * The derivatives of the bias potential imposed on the system.
         * They are evaluated by using the lookup table.
         */
		std::vector<double> derivatives_;

        //! The order of the basis polynomials
        std::vector<unsigned int> polyords_;

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

		//! Updates the bias projection of the PMF.
        /*!
         * \param cvs List of collective variables.
         * \param beta Temperature equivalent.
         */
		void UpdateBias(const CVList& cvs, const double);

		//! Computes the bias force.
        /*!
         * \param cvs List of collective variables.
         */
		void CalcBiasForce(const CVList& cvs);

		//! Prints the current bias to a defined file from the JSON.
        /*!
         * \param cvs List of collective variables.
         * \param beta Scale parameter.
         */
		void PrintBias(const CVList& cvs, const double beta);

        //! Initializes the look up tables for polynomials
        /*!
         * \param cvs List of collective variables.
         */
        void BasisInit(const CVList& cvs);

		//! Output stream for basis projection data.
		std::ofstream basisout_;

        //! Output stream for coefficients (for reading purposes)
        std::ofstream coeffout_;

        //! The option to name both the basis and coefficient files will be given
        //! Basis filename 
        std::string bnme_;

        //! Coefficient filename
        std::string cnme_;

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
		Basis(const MPI_Comm& world,
			  const MPI_Comm& comm,
              Histogram<int> *hist,
			  const std::vector<unsigned int>& polyord,
              const std::vector<double>& restraint,
              const std::vector<double>& boundUp,
              const std::vector<double>& boundLow,
              unsigned int cyclefreq,
			  unsigned int frequency,
              const std::string bnme,
              const std::string cnme,
              const double temperature,
              const double tol,
              const double weight,
              bool converge) :
		Method(frequency, world, comm), hist_(hist),
        coeff_(), unbias_(), coeff_arr_(), LUT_(), derivatives_(), polyords_(polyord),
        restraint_(restraint), boundUp_(boundUp), boundLow_(boundLow),
        cyclefreq_(cyclefreq), mpiid_(0), weight_(weight),
        temperature_(temperature), tol_(tol),
        converge_exit_(converge), bnme_(bnme), cnme_(cnme), iteration_(0)
		{
		}

		//! Pre-simulation hook.
        /*!
         * \param snapshot Simulation snapshot.
         * \param cvs List of CVs.
         */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
        /*!
         * \param snapshot Simulation snapshot.
         * \param cvs List of CVs.
         */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
        /*!
         * \param snapshot Simulation snapshot.
         * \param cvs List of CVs.
         */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

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
		static Basis* Construct(const Json::Value& json, 
		                        const MPI_Comm& world,
		                        const MPI_Comm& comm,
					            const std::string& path);

        //! \copydoc Serializable::Serialize()
        /*!
         * \warning Serialization is not implemented yet!
         */
		void Serialize(Json::Value& json) const override;

        //! Destructor.
        ~Basis()
        {
            delete hist_;
        }
	};
}
			
