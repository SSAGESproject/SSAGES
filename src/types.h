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
