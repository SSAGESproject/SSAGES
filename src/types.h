#pragma once 
#include <Eigen/Dense>

namespace SSAGES
{    
	//! Three-dimensional vector.
	using Vector3 = Eigen::Vector3d;

	//! Three-dimensional boolean.
	using Bool3 = Eigen::Matrix<bool, 3, 1>;

	//! 3x3 matrix. 
	using Matrix3 = Eigen::Matrix3d;

	//! List of integers.
	using Label = std::vector<int>;
}