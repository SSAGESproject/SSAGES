#pragma once

#include "Grid.h"

#include "../JSON/Serializable.h"
#include <vector>

namespace SSAGES
{
	// Interface for a collective variable.
	class Grid1D : public Grid
	{
	private:

		double _lower;
		double _upper;
		bool _periodic;
		int _num_points;
		double _spacing;
		std::vector<float> _values;

	public:

		Grid1D(double lower, double upper,
			bool periodic, int num_points) : 
		_lower(lower), _upper(upper), _periodic(periodic), _num_points(num_points) 
		{

			//Generate Grid
			_spacing = (_upper - _lower)/double(_num_points - 1);

			for(int i=0;i<num_points; i++)
				_values.push_back(0);

		}
	
		// Return the nearest index for given values.
		// Out of bounds will return edge vertex
		std::vector<int> GetIndices(const std::vector<float> &val) const override
		{
			std::vector<int> vertices;

			for(size_t i = 0; i < val.size(); i++)
			{
				int vertex = int(val[i] - _lower + (_spacing/2.0));
				if(vertex < 0) // out of bounds
					vertex = 0;
				else if(vertex > _num_points -1) // out of bounds
				{
					if(_periodic)
						vertex = 0;
					else
						vertex = _num_points -1;
				}
				vertices.push_back(vertex);
			}

			return vertices;
		}

		// Return a pointer to the grid value allowing user to use it
		// or modify it.
		float GetValue(const std::vector<int>& indices) const override
		{
			return _values[indices[0]];
		}

		void SetValue(const std::vector<int>& indices, float value) override
		{
			_values[indices[0]] = value;
		}

		//Currently only for debugging
		void PrintGrid()
		{
			for(size_t i = 0; i<_values.size(); i++)
			{
				std::cout << _values[i]<<" ";
			}
			std::cout<<std::endl;
		}

	};	
}