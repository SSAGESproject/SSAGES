#pragma once

#include "Grid.h"

#include "../JSON/Serializable.h"
#include <vector>

namespace SSAGES
{
	//! 3D Grid.
	class Grid3D : public Grid
	{
	private:

		//!< 3D Vector containing the Grid values.
		std::vector<std::vector<std::vector<float>>> _values;
		std::vector<std::vector<std::vector<std::array<double, 3>>>> _derivs;

	public:

		//! Constructor.
		/*!
		 * \param lower List of values for the lower edges of the Grid.
		 * \param upper List of values for the upper edges of the Grid.
		 * \param periodic List of periodicity values.
		 * \param num_points List of how many Grid points in each dimension.
		 */
		Grid3D(std::vector<double> lower, std::vector<double> upper,
			std::vector<bool> periodic, std::vector<int> num_points)
		{
			_NDim = 3;
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
			  _values[i].resize(_num_points[1]);
			  _derivs[i].resize(_num_points[1]);
			}
						

			for(size_t i = 0; i < _values.size(); i++){
			  for(size_t j = 0; j < _values[i].size(); j++){
			    _values[i][j].resize(_num_points[2],0);
			    _derivs[i][j].resize(_num_points[2],{0,0,0});
			  }
			}
		}

		//! Get value at a given Grid point.
		/*!
		 * \param indices List of indices specifying the Grid point.
		 * \return Value at this grid point.
		 */
		float GetValue(const std::vector<int>& indices) const override
		{
			return _values[indices[0]][indices[1]][indices[2]];
		}
		float GetDeriv(const std::vector<int>& indices, int dim) const override
		{
		  return _derivs[indices[0]][indices[1]][indices[2]][dim];
		}

		float InterpolateValue(const std::vector<float> &val){
		  std::vector<int> vertices;
		  std::vector<float> gridpos;
		  std::vector<float> gridval;
		  int ii[3];
		  float ival;

		  for(ii[2] = 0; ii[2] <= 1; ii[2]++){
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
		  }

		  //now, do 3d interpolation
		  //order of vertices is 0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 1,0,1; 0,1,1; 1,1,1;
		  ival  = 0;
		  ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1])*(val[2] - gridpos[0][2])*gridval[7];
		  ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][1])*(val[2] - gridpos[1][2])*gridval[6];
		  ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1])*(val[2] - gridpos[2][2])*gridval[5];
		  ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1])*(val[2] - gridpos[3][2])*gridval[4];
		  ival += -(val[0] - gridpos[4][0])*(val[1] - gridpos[4][1])*(val[2] - gridpos[4][2])*gridval[3];
		  ival +=  (val[0] - gridpos[5][0])*(val[1] - gridpos[5][1])*(val[2] - gridpos[5][2])*gridval[2];
		  ival +=  (val[0] - gridpos[6][0])*(val[1] - gridpos[6][1])*(val[2] - gridpos[6][2])*gridval[1];
		  ival += -(val[0] - gridpos[7][0])*(val[1] - gridpos[7][1])*(val[2] - gridpos[7][2])*gridval[0];

		  
		  ival /= _grid->spacing[0]*_grid->spacing[1]_grid->spacing[2;
			
		  return ival;
		}

		    float InterpolateDeriv(const std::vector<float> &val, int dim){
		  std::vector<int> vertices;
		  std::vector<float> gridpos;
		  std::vector<float> gridval;
		  int ii[3];
		  float ival;

		  for(ii[2] = 0; ii[2] <= 1; ii[2]++){
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
		  }

		  //now, do 3d interpolation
		  //order of vertices is 0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 1,0,1; 0,1,1; 1,1,1;
		  ival  = 0;
		  ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1])*(val[2] - gridpos[0][2])*gridval[7];
		  ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][1])*(val[2] - gridpos[1][2])*gridval[6];
		  ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1])*(val[2] - gridpos[2][2])*gridval[5];
		  ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1])*(val[2] - gridpos[3][2])*gridval[4];
		  ival += -(val[0] - gridpos[4][0])*(val[1] - gridpos[4][1])*(val[2] - gridpos[4][2])*gridval[3];
		  ival +=  (val[0] - gridpos[5][0])*(val[1] - gridpos[5][1])*(val[2] - gridpos[5][2])*gridval[2];
		  ival +=  (val[0] - gridpos[6][0])*(val[1] - gridpos[6][1])*(val[2] - gridpos[6][2])*gridval[1];
		  ival += -(val[0] - gridpos[7][0])*(val[1] - gridpos[7][1])*(val[2] - gridpos[7][2])*gridval[0];

		  
		  ival /= _grid->spacing[0]*_grid->spacing[1]_grid->spacing[2;
			
		  return ival;
		}
		    
		//! Set value at a given Grid point.
		/*!
		 * \param indices List of indices specifying the Grid point.
		 * \param value New value for the specified Grid point.
		 */
		void SetValue(const std::vector<int>& indices, float value) override
		{
			_values[indices[0]][indices[1]][indices[2]] = value;
		}
		void SetDeriv(const std::vector<int>& indices, float value, int dim) override
		{
			_derivs[indices[0]][indices[1]][indices[2]][dim] = value;
		}

		//! Write the Grid to console output.
		/*!
		 * Currently only for debugging.
		 */
		void PrintGrid() const override
		{
			for(size_t i = 0; i<_values.size(); i++)
			{
				for(size_t j = 0; j< _values[i].size(); j++)
				{
					for(size_t k = 0; k < _values[i][j].size(); k++)
						std::cout << _values[i][j][k]<<" ";
					std::cout<<std::endl;
				}
				std::cout<<std::endl;
			}
			std::cout<<std::endl;
		}
	};	
}
