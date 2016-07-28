#pragma once

#include "Grid.h"

#include "../JSON/Serializable.h"
#include <vector>
#include <iterator>

namespace SSAGES
{
	//! 1D Grid.
	class Grid1D : public Grid
	{
	private:

	public:

		//! Constructor
		/*!
		 * \param lower List of lower edge values.
		 * \param upper List of upper edge values.
		 * \param periodic List of Grid periodicity values.
		 * \param num_points List of numbers of Grid points.
		 */
		Grid1D(std::vector<double> lower, std::vector<double> upper,
			std::vector<bool> periodic, std::vector<int> num_points)
		{
			_NDim = 1;
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
			_flatvector.resize(_num_points[0]);
			for(int i = 0; i < _num_points[0]; i++)
			{
				_flatvector[i].first = 0.0;
				_flatvector[i].second.push_back(_lower[i] + _spacing[0]*i);
			}
		}

		~Grid1D(){}
	};	
}