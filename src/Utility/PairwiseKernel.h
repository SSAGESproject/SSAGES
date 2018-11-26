/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
 *                Michael Quevillon <mquevill@nd.edu>
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

#include <cmath>
#include "json/json.h"

namespace SSAGES
{
	//! Pairwise kernel base class.
	/*!
	 *
	 * Pairwise functions are useful for a variety of collective variables.
	 * Several classes of these functions exist.
	 * 
	 * Gaussian functions are useful for collective variables that focus on
	 * a certain value. Unlike a delta function, they interpolate provide 
	 * a continuous function and gradient for biasing. 
	 * 
	 * Switching functions are useful for collective variables that are discrete 
	 * in nature. They interpolate the discrete values and provide a continuous 
	 * gradient for biasing.
	 * 
	 * \ingroup Core
	 */
	class PairwiseKernel
	{
	public:

		//! Evaluate the pairwise kernel function.
		/*!
		 * \param rij distance between two atoms. 
		 * \param df Reference to variable which will store the gradient.
		 * 
		 * \return value of pairwise kernel function. 
		 */
		virtual double Evaluate(double rij, double& df) const = 0;

		//! Build PairwiseKernel from JSON value. 
		/*!
		 * \param json JSON value node.
		 * \param path Path for JSON path specification.
		 * 
		 * \return Pointer to new PairwiseKernel.
		 */
		static PairwiseKernel* Build(const Json::Value& json, const std::string& path);
		
		//! Destructor
		virtual ~PairwiseKernel() {}
	};

	//! Gaussian Function 
	/*!
	 *
	 * A standard Gaussian function (also called a "bell curve").
	 * 
	 */
	class GaussianPK : public PairwiseKernel
	{
	private: 
		double mu_; //!< Center of Gaussian.
		double sigma_; //!< Width of Gaussian.

	public:
		//! Constructor.
		/*!
		 * \param mu Center of Gaussian.
		 * \param sigma Width of Gaussian.
		 *
		 * Construct a GaussianPK.
		 *
		 */
		GaussianPK(double mu, double sigma) : 
		mu_(mu), sigma_(sigma) {}
		
		double Evaluate(double rij, double& df) const
		{
			const auto dx = (rij - mu_)/sigma_;
			const auto f = exp( - dx*dx/2.);
			const auto pre = - dx/sigma_;
			
			df = pre * f;
			return 	f;
		}

		//! Build GaussianPK from JSON value. 
		/*!
		 * \param json JSON value node.
		 * \param path Path for JSON path specification.
		 * 
		 * \return Pointer to new GaussianPK.
		 */
		static GaussianPK* Build(const Json::Value& json, const std::string& path);
	};

	//! Rational Switching Function 
	/*!
	 *
	 * A switching function of a rational form.
	 * 
	 */
	class RationalSwitchPK : public PairwiseKernel
	{
	private: 
		double d0_; //!< Minimum linear shift value. 
		double r0_; //!< Cutoff distance. 
		int n_;     //!< Exponent of numerator in switching function (controls stiffness).
		int m_;     //!< Exponent of denominator in switching function (controls stiffness). 
	
	public:
		//! Constructor.
		/*!
		 * \param d0 Minimum linear shift value.
		 * \param r0 Cutoff distance.
		 * \param m An exponent of the switching function.
		 * \param n An exponent of the switching function.
		 *
		 * Construct a GaussianPK.
		 *
		 */
		RationalSwitchPK(double d0, double r0, int n, int m) : 
		d0_(d0), r0_(r0), n_(n), m_(m) {}

		double Evaluate(double rij, double& df) const override
		{
			const auto xarg = (rij - d0_)/r0_;
			const auto xn = std::pow(xarg, n_);
			const auto xm = std::pow(xarg, m_);
			const auto f = (1.-xn)/(1.-xm);
			
			df = f/(d0_-rij)*(n_*xn/(1.-xn)+m_*xm/(xm-1.));
			return 	f;
		}

		//! Build RationalSwitchPK from JSON value. 
		/*!
		 * \param json JSON value node.
		 * \param path Path for JSON path specification.
		 * 
		 * \return Pointer to new RationalSwitchPK.
		 */
		static RationalSwitchPK* Build(const Json::Value& json, const std::string& path);
	};
}