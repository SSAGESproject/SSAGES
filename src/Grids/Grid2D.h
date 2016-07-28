#pragma once

#include "Grid.h"

#include "../JSON/Serializable.h"
#include <vector>

namespace SSAGES
{
	//! 2D Grid.
	class Grid2D : public Grid
	{
	private:

	public:

		//! Constuctor
		/*!
		 * \param lower List of values for the lower edges of the Grid.
		 * \param upper List of values for the upper edges of the Grid.
		 * \param periodic List of periodicity values.
		 * \param num_points List of how many Grid points in each direction.
		 */
		Grid2D(std::vector<double> lower, std::vector<double> upper,
			std::vector<bool> periodic, std::vector<int> num_points)
		{
			_NDim = 2;
			for(size_t i =0; i <lower.size(); i++)
			{
				_lower.push_back(lower[i]);
				_upper.push_back(upper[i]);
				_periodic.push_back(periodic[i]);
				_num_points.push_back(num_points[i]);
				_spacing.push_back(0.0);
			}
			//Generate Grid
			for(size_t i = 0; i < _spacing.size(); i++)
				_spacing[i] = (_upper[i] - _lower[i])/double(_num_points[i] - 1);

			// Construct flat vector 
			_flatvector.resize(_num_points[0]*_num_points[1]);
			for(int i = 0; i < _num_points[0]; i++)
			{
				for(int j = 0; j < _num_points[1]; j++)
				{
					int flat = FlattenIndices({i,j},_num_points);
					_flatvector[flat].first = 0;
					_flatvector[flat].second.push_back(_lower[i] + _spacing[0]*i);
					_flatvector[flat].second.push_back(_lower[j] + _spacing[1]*j);
				}
			}
		}

		~Grid2D(){}
	};	
}