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
					int flat = FlattenIndices({i,j}, _num_points);
					_flatvector[flat].first = std::vector<double>(3, 0.0);
					_flatvector[flat].second.push_back(_lower[0] + _spacing[0]*i);
					_flatvector[flat].second.push_back(_lower[1] + _spacing[1]*j);
				}
			}
		}

		std::vector<std::vector<int>> GetVoxel(const std::vector<double> &val) const override
		{
			std::vector<int> vertices;
			std::vector<std::vector<int>> voxel;

			for(size_t i = 0; i < val.size(); i++)
			{
				int vertex = 0;

				vertex = int((val[i] - _lower[i])/_spacing[i]);

				if(_periodic[i])
				{
					vertex = vertex % _num_points[i];
					if(vertex<0)
						vertex += _num_points[i];
				}

				if(vertex < 0 || vertex >=_num_points[i]) // out of bounds
					throw std::out_of_range("Voxel not in grid!");

				vertices.push_back(vertex);
			}

			//order of vertices is 0,0; 1,0; 0,1; 1,1
			int oldj = vertices[0];
			for(size_t i = 0; i < 2; i++)
			{
				vertices[1] += i;
				vertices[0] = oldj;
				for(size_t j =0; j < 2; j++)
				{
					vertices[0] += j;
					voxel.push_back(vertices);
				}
			}

			return voxel;
		}

		double InterpolateValue(const std::vector<double> &val) const override
		{
			std::vector<std::vector<int>> voxel = GetVoxel(val);
			std::vector<std::vector<double>> gridpos;
			std::vector<double> gridval;

		  	for(size_t v = 0; v < voxel.size(); v++)
		  	{
		  		gridpos.push_back(GetLocation(voxel[v]));
		  		gridval.push_back(GetValue(voxel[v]));
		  	}

			//now, do 2d interpolation
			//order of vertices is 0,0; 1,0; 0,1; 1,1
			//note that middle two are negated so arguments are positive, since val < gridpos in
			// one of two directions.
			double ival  = 0;
			ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1]) * gridval[3];
			ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][0]) * gridval[2];
			ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1]) * gridval[1];
			ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1]) * gridval[0];


			ival /= _spacing[0]*_spacing[1];

			return ival;
		}

		double InterpolateDeriv(const std::vector<double> &val, int dim) const override
		{
			std::vector<std::vector<int>> voxel = GetVoxel(val);
			std::vector<std::vector<double>> gridpos;
			std::vector<std::vector<double>> gridval;

		  	for(size_t v = 0; v < voxel.size(); v++)
		  	{
		  		gridpos.push_back(GetLocation(voxel[v]));
		  		gridval.push_back(GetExtra(voxel[v]));
		  	}

			//now, do 2d interpolation
			//order of vertices is 0,0; 1,0; 0,1; 1,1
			//note that middle two are negated so arguments are positive, since val < gridpos in
			// one of two directions.
			double ival  = 0;
			ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1]) * gridval[3][dim];
			ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][0]) * gridval[2][dim];
			ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1]) * gridval[1][dim];
			ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1]) * gridval[0][dim];

			ival /= _spacing[0]*_spacing[1];

			return ival;
		}

		~Grid2D(){}
	};	
}