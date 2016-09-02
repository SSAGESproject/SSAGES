/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
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

	public:

		//! Constructor.
		/*!
		 * \param lower List of values for the lower edges of the Grid.
		 * \param upper List of values for the upper edges of the Grid.
		 * \param num_points List of how many Grid points in each dimension.
		 */
		Grid3D(std::vector<double> lower, std::vector<double> upper, std::vector<int> num_points)
		{
			_NDim = 3;
			for(size_t i =0; i <lower.size(); i++)
			{
				_lower.push_back(lower[i]);
				_upper.push_back(upper[i]);
				_periodic.push_back(false);
				_num_points.push_back(num_points[i]);
				_spacing.push_back(0.0);
			}
			//Generate Grid
			for(size_t i = 0; i < _spacing.size(); i++)
				_spacing[i] = (_upper[i] - _lower[i])/double(_num_points[i] - 1);

			// Construct flat vector 
			_flatvector.resize(_num_points[0]*_num_points[1]*_num_points[2]);
			for(int i = 0; i < _num_points[0]; i++)
			{
				for(int j = 0; j < _num_points[1]; j++)
				{
					for(int k = 0; k < _num_points[2]; k++)
					{
						int flat = FlattenIndices({i,j,k},_num_points);
						_flatvector[flat].first = std::vector<double>(4, 0.0);
						_flatvector[flat].second.push_back(_lower[0] + _spacing[0]*i);
						_flatvector[flat].second.push_back(_lower[1] + _spacing[1]*j);
						_flatvector[flat].second.push_back(_lower[2] + _spacing[2]*k);
					}
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

			//order of vertices is 0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 1,0,1; 0,1,1; 1,1,1;
			int oldj = vertices[1];
			int oldk = vertices[0];
			for(size_t i = 0; i < 2; i++)
			{
				vertices[2] += i;
				vertices[1] = oldj;
				for(size_t j =0; j < 2; j++)
				{
					vertices[1] += j;
					vertices[0] = oldk;
					for(size_t k = 0; k < 2; k++)
					{
						vertices[0] += k;
						voxel.push_back(vertices);
					}
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

			//now, do 3d interpolation
			//order of vertices is 0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 1,0,1; 0,1,1; 1,1,1;
			double ival  = 0;
			ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1])*(val[2] - gridpos[0][2])*gridval[7];
			ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][1])*(val[2] - gridpos[1][2])*gridval[6];
			ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1])*(val[2] - gridpos[2][2])*gridval[5];
			ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1])*(val[2] - gridpos[3][2])*gridval[4];
			ival += -(val[0] - gridpos[4][0])*(val[1] - gridpos[4][1])*(val[2] - gridpos[4][2])*gridval[3];
			ival +=  (val[0] - gridpos[5][0])*(val[1] - gridpos[5][1])*(val[2] - gridpos[5][2])*gridval[2];
			ival +=  (val[0] - gridpos[6][0])*(val[1] - gridpos[6][1])*(val[2] - gridpos[6][2])*gridval[1];
			ival += -(val[0] - gridpos[7][0])*(val[1] - gridpos[7][1])*(val[2] - gridpos[7][2])*gridval[0];


			ival /= _spacing[0]*_spacing[1]*_spacing[2];

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

			//now, do 3d interpolation
			//order of vertices is 0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 1,0,1; 0,1,1; 1,1,1;
			double ival  = 0;
			ival +=  (val[0] - gridpos[0][0])*(val[1] - gridpos[0][1])*(val[2] - gridpos[0][2])*gridval[7][dim];
			ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][1])*(val[2] - gridpos[1][2])*gridval[6][dim];
			ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1])*(val[2] - gridpos[2][2])*gridval[5][dim];
			ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1])*(val[2] - gridpos[3][2])*gridval[4][dim];
			ival += -(val[0] - gridpos[4][0])*(val[1] - gridpos[4][1])*(val[2] - gridpos[4][2])*gridval[3][dim];
			ival +=  (val[0] - gridpos[5][0])*(val[1] - gridpos[5][1])*(val[2] - gridpos[5][2])*gridval[2][dim];
			ival +=  (val[0] - gridpos[6][0])*(val[1] - gridpos[6][1])*(val[2] - gridpos[6][2])*gridval[1][dim];
			ival += -(val[0] - gridpos[7][0])*(val[1] - gridpos[7][1])*(val[2] - gridpos[7][2])*gridval[0][dim];


			ival /= _spacing[0]*_spacing[1]*_spacing[2];

			return ival;
		}

		virtual ~Grid3D(){}
	};	
}