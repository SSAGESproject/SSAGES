#pragma once

#include "Grid.h"

#include "../JSON/Serializable.h"
#include <vector>

namespace SSAGES
{
	//! 1D Grid.
	class Grid1D : public Grid
	{
	private:

		std::vector<double> _values; //!< Grid values
		std::vector<std::array<double, 1>> _derivs; //!< Grid derivatives
		
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

			_values.resize(_num_points[0],0);
			_derivs.resize(_num_points[0],{0});
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
		
		float GetDeriv(const std::vector<int>& indices, int dim) const override
		{
		  return _derivs[indices[0]][dim];
		}

		float InterpolateValue(const std::vector<float> &val){
		  std::vector<int> vertices;
		  std::vector<float> gridpos;
		  std::vector<float> gridval;
		  int ii;
		  float ival;
		  
		  for(ii = 0; ii <= 1; ii++){
		    for(size_t i = 0; i < val.size(); i++)
		      {
			//always round down, allows us get the vertices appropriately
			int vertex = int((val[i] - _lower[i])/_spacing[i]) + ii;
			if(vertex < 0) // out of bounds
			  vertex = 0;
			else if(vertex > _num_points[i] -1) // out of bounds
			{
			  if(_periodic[i])
			    vertex = 0;
			  else
			    vertex = _num_points[i] -1;
			}
		      vertices.push_back(vertex);
		    }
		    gridpos=push_back(GetLocation(vertices));
		    gridval=push_back(GetValue(vertices));
		    vertices.clear();
		  }

		  //now, do 1d interpolation
		  ival = ((val[i]-gridpos[0])*gridval[1] +
			  (gridpos[1]-val[i])*gridval[0]) /
		          (gridpos[1]-gridpos[0]);

		  return ival;
		}

		float InterpolateDeriv(const std::vector<float> &val, int dim){
		  std::vector<int> vertices;
		  std::vector<float> gridpos;
		  std::vector<float> gridval;
		  int ii;
		  float ival;
		  
		  for(ii = 0; ii <= 1; ii++){
		    for(size_t i = 0; i < val.size(); i++)
		      {
			//always round down, allows us get the vertices appropriately
			int vertex = int((val[i] - _lower[i])/_spacing[i]) + ii;
			if(vertex < 0) // out of bounds
			  vertex = 0;
			else if(vertex > _num_points[i] -1) // out of bounds
			{
			  if(_periodic[i])
			    vertex = 0;
			  else
			    vertex = _num_points[i] -1;
			}
		      vertices.push_back(vertex);
		    }
		    gridpos=push_back(GetLocation(vertices));
		    gridval=push_back(GetDeriv(vertices,dim));
		    vertices.clear();
		  }

		  //now, do 1d interpolation
		  ival = ((val[i]-gridpos[0])*gridval[1] +
			  (gridpos[1]-val[i])*gridval[0]) /
		          (gridpos[1]-gridpos[0]);

		  return ival;
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
		void SetDeriv(const std::vector<int>& indices, float value, int dim) override
		{
			_derivs[indices[0]][dim] = value;
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

	};	
}
