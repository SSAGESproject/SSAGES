/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
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
#include <Eigen/Dense>
#include <vector>

namespace SSAGES
{    
	// Forward declare necessary types..
	class CollectiveVariable;

	//! Three-dimensional vector.
	using Vector3 = Eigen::Vector3d;

	//! Three-dimensional boolean.
	using Bool3 = Eigen::Matrix<bool, 3, 1>;

	//! Three-dimensional integer vector.
	using Integer3 = Eigen::Vector3i;

	//! 3x3 matrix. 
	using Matrix3 = Eigen::Matrix3d;

	//! List of integers.
	using Label = std::vector<int>;

	//! List of Collective Variables.
	using CVList = std::vector<CollectiveVariable*>;
}
