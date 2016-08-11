#pragma once 
#include <Eigen/Dense>

namespace SSAGES
{    
	// Forward declare necessary types..
	class CollectiveVariable;

	//! Three-dimensional vector.
	using Vector3 = Eigen::Vector3d;

	//! 3x3 matrix. 
	using Matrix3 = Eigen::Matrix3d;

	//! List of integers.
	using Label = std::vector<int>;

	//! List of Collective Variables.
	using CVList = std::vector<CollectiveVariable*>;
}