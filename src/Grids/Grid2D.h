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

		//! Vector of vector containing values.
		std::vector<std::vector<double>> _values;

	public:

		//using const_iterator = Grid2D::const_iterator;

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

			_values.resize(_num_points[0]);

			for(size_t i = 0; i < _values.size(); i++)
				_values[i].resize(_num_points[1],0);

			// Construct flat vector 
			_flatvector.resize(_num_points[0]*_num_points[1]);
			for(int i = 0; i < _num_points[0]; i++)
			{
				for(int j = 0; j < _num_points[1]; j++)
				{
					auto loc = GetLocation({i,j});
					_flatvector[i*_num_points[1] + j].first[0] = loc[0];
					_flatvector[i*_num_points[1] + j].first[1] = loc[1];
					_flatvector[i*_num_points[1] + j].second = _values[i][j];
				}
			}
		}
	
		//! Get the value at a given Grid point.
		/*!
		 * \param indices List of indices specifying the Grid point.
		 * \return Value at the specified Grid point.
		 */
		float GetValue(const std::vector<int>& indices) const override
		{
			return _values[indices[0]][indices[1]];
		}

		//! Set the value at a given Grid point.
		/*!
		 * \param indices List of indices specifying the Grid point.
		 * \param value New value for the specified Grid point.
		 */
		void SetValue(const std::vector<int>& indices, float value) override
		{
			_values[indices[0]][indices[1]] = value;
		}

		//! Write grid to console
		/*!
		 * Currently only for debugging.
		 */
		void PrintGrid() const override
		{
			for(size_t i = 0; i<_values.size(); i++)
			{
				for(size_t j = 0; j< _values[i].size(); j++)
					std::cout << _values[i][j]<<" ";
				std::cout<<std::endl;
			}
			std::cout<<std::endl;
		}

		//Grid2D::const_iterator begin() const { _flatvector.begin();}
		//Grid2D::const_iterator end() const {_flatvector.end();}

	};	
}