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

		std::vector<double> _values; //!< Grid values

	public:

		//using const_iterator = Grid1D::const_iterator;

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

			_values.resize(_num_points[0]);


			// Construct flat vector 
			_flatvector.resize(_num_points[0]);
			for(int i = 0; i < _num_points[0]; i++)
			{
				auto loc = GetLocation({i});
				_flatvector[i].first[0] = loc[0];
				_flatvector[i].second = _values[i];
			}
		}

		//! Get value at one Grid point.
		/*!
		 * \param indices List of indices. Must be at least of size 1.
		 * \return Value at the grid point specified by the first element of indices.
		 */
		float GetValue(const std::vector<int>& indices) const override
		{
			return _values[indices[0]];
		}

		//! Set value at one grid point.
		/*!
		 * \param indices List of indices. Must be at least of size 1.
		 * \param value New value for the grid point specified by the first element in indices.
		 */
		void SetValue(const std::vector<int>& indices, float value) override
		{
			_values[indices[0]] = value;
		}

		//! Write Grid to console
		/*!
		 * Currently only for debugging
		 */
		void PrintGrid() const override
		{
			for(size_t i = 0; i<_values.size(); i++)
			{
				std::cout << _values[i]<<" ";
			}
			std::cout<<std::endl;
		}

		//Grid1D::const_iterator begin() const { _flatvector.begin();}
		//Grid1D::const_iterator end() const {_flatvector.end();}

	};	
}