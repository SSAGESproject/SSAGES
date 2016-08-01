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
		std::vector<std::vector<std::array<double, 2>>> _derivs;
		
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

			_values.resize(_num_points[0]);
			_derivs.resize(_num_points[0]);

			for(size_t i = 0; i < _values.size(); i++){
			  _values[i].resize(_num_points[1],0);
			  _derivs[i].resize(_num_points[1],{0,0});
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
		float GetDeriv(const std::vector<int>& indices, int dim) const override
		{
		  return _derivs[indices[0]][indices[1]][dim];
		}

		float InterpolateValue(const std::vector<float> &val){
		  std::vector<int> vertices;
		  std::vector<float> gridpos;
		  std::vector<float> gridval;
		  int ii[2];
		  float ival;

		  for(ii[1] = 0; ii[1] <= 1; ii[1]++){
		    for(ii[0] = 0; ii[0] <= 1; ii[0]++){
		      for(size_t i = 0; i < val.size(); i++)
			{
			  //always round down, allows us get the vertices appropriately
			  int vertex = int((val[i] - _lower[i])/_spacing[i]) + ii[i];
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
		  }

		  //now, do 2d interpolation
		  //order of vertices is 0,0; 1,0; 0,1; 1,1
		  //note that middle two are negated so arguments are positive, since val < gridpos in
		  // one of two directions.
		  ival  = 0;
		  ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1]) * gridval[3];
		  ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][0]) * gridval[2];
		  ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1]) * gridval[1];
		  ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1]) * gridval[0];

		  
		  ival /= _grid->spacing[0]*_grid->spacing[1];
			
		  return ival;
		}

		float InterpolateDeriv(const std::vector<float> &val, int dim){
		  std::vector<int> vertices;
		  std::vector<float> gridpos;
		  std::vector<float> gridval;
		  int ii[2];
		  float ival;
		  
		  for(ii[1] = 0; ii[1] <= 1; ii[1]++){
		    for(ii[0] = 0; ii[0] <= 1; ii[0]++){
		      for(size_t i = 0; i < val.size(); i++)
			{
			  //always round down, allows us get the vertices appropriately
			  int vertex = int((val[i] - _lower[i])/_spacing[i]) + ii[i];
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
		  }

		  //now, do 2d interpolation
		  //order of vertices is 0,0; 1,0; 0,1; 1,1
		  //note that middle two are negated so arguments are positive, since val < gridpos in
		  // one of two directions.
		  ival  = 0;
		  ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1]) * gridval[3];
		  ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][0]) * gridval[2];
		  ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1]) * gridval[1];
		  ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1]) * gridval[0];

		  
		  ival /= _grid->spacing[0]*_grid->spacing[1];
		  return ival;
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
		void SetDeriv(const std::vector<int>& indices, float value, int dim) override
		{
			_derivs[indices[0]][indices[1]][dim] = value;
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

	};	
}
