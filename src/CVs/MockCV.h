/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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

#include "CollectiveVariable.h"

#include <array>
#include <cmath>

namespace SSAGES
{
	//! Mock collective variable for testing purposes.
	class MockCV : public CollectiveVariable
	{
	private:
		//! User defined gradient vector.
		Vector3 _user_defined_grad;

	public:

		//! Constructor.
		/*!
		 * \param value Value the Mock CV will return.
		 * \param grad Gradient vector the Mock CV will return.
		 * \param upper Value for the upper bound of the CV.
		 * \param lower Value for the lower bound of the CV.
		 *
		 * Construct a mock CV
		 */
		MockCV(double value, const Vector3& grad, double upper, double lower) :
		_user_defined_grad(grad)
		{
			_val = value;
			_bounds = {{upper, lower}};
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient
			auto n = snapshot.GetPositions().size();
			_grad.resize(n, _user_defined_grad);
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{

		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		void Serialize(Json::Value& json) const override
		{

		}

	};
}
