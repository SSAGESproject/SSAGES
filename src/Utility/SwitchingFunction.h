/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
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

#include <cmath>
#include "json/json.h"

namespace SSAGES
{
	//! Switching function base class. 
	/*!
	 *
	 * Switching functions are useful for collective variables that are discrete 
	 * in nature. They interpolate the discrete values and provide a continuous 
	 * gradient for biasing. 
	 * 
	 * \ingroup Core
	 */
	class SwitchingFunction
	{
	public:

		//! Evaluate the switching function.
		/*!
		 * \param rij distance between two atoms. 
		 * \param df Reference to variable which will store the gradient.
		 * 
		 * \return value of switching function. 
		 */
		virtual double Evaluate(double rij, double& df) const = 0;

		//! Build SwitchingFunction from JSON value. 
		/*!
		 * \param json JSON value node. 
		 * 
		 * \return Pointer to new SwitchingFunction.
		 */
		static SwitchingFunction* Build(const Json::Value& json);
	};

	class RationalSF : public SwitchingFunction
	{
	private: 
		double d0_; //!< Minimum linear shift value. 
		double r0_; //!< Cutoff distance. 
		int n_, m_; //!< Exponents of the switching function which control the stiffness. 
	
	public:
		RationalSF(double d0, double r0, int n, int m) : 
		d0_(d0), r0_(r0), n_(n), m_(m) {}

		//! Evaluate the switching function.
		/*!
		 * \param rij distance between two atoms. 
		 * \param df Reference to variable which will store the gradient.
		 * 
		 * \return value of switching function. 
		 */
		double Evaluate(double rij, double& df) const override
		{
			const auto xarg = (rij - d0_)/r0_;
			const auto xn = std::pow(xarg, n_);
			const auto xm = std::pow(xarg, m_);
			const auto f = (1.-xn)/(1.-xm);
			
			df = f/(d0_-rij)*(n_*xn/(1.-xn)+m_*xm/(xm-1.));
			return 	f;
		}

		//! Build SwitchingFunction from JSON value. 
		/*!
		 * \param json JSON value node. 
		 * 
		 * \return Pointer to new SwitchingFunction.
		 */
		static SwitchingFunction* Build(const Json::Value& json)
		{
			return new RationalSF(
						json["d0"].asDouble(), 
						json["r0"].asDouble(), 
						json["n"].asInt(), 
						json["m"].asInt()
					);
		}
	};
}