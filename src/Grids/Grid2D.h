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
	//! 2D Grid.
	class Grid2D : public Grid
	{
	private:

	public:

		//! Constuctor
		/*!
		 * \param lower List of values for the lower edges of the Grid.
		 * \param upper List of values for the upper edges of the Grid.
		 * \param num_points List of how many Grid points in each direction.
		 */
		Grid2D(std::vector<double> lower, std::vector<double> upper, std::vector<int> num_points)
		{
			NDim_ = 2;
			for(size_t i =0; i <lower.size(); i++)
			{
				lower_.push_back(lower[i]);
				upper_.push_back(upper[i]);
				periodic_.push_back(false);
				num_points_.push_back(num_points[i]);
				spacing_.push_back(0.0);
			}
			//Generate Grid
			for(size_t i = 0; i < spacing_.size(); i++)
				spacing_[i] = (upper_[i] - lower_[i])/double(num_points_[i] - 1);

			// Construct flat vector 
			flatvector_.resize(num_points_[0]*num_points_[1]);
			for(int i = 0; i < num_points_[0]; i++)
			{
				for(int j = 0; j < num_points_[1]; j++)
				{
					int flat = FlattenIndices({i,j}, num_points_);
					flatvector_[flat].first = std::vector<double>(3, 0.0);
					flatvector_[flat].second.push_back(lower_[0] + spacing_[0]*i);
					flatvector_[flat].second.push_back(lower_[1] + spacing_[1]*j);
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

				vertex = int((val[i] - lower_[i])/spacing_[i]);

				if(periodic_[i])
				{
					vertex = vertex % num_points_[i];
					if(vertex<0)
						vertex += num_points_[i];
				}

				if(vertex < 0 || vertex >=num_points_[i]) // out of bounds
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
			ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][1]) * gridval[2];
			ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1]) * gridval[1];
			ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1]) * gridval[0];


			ival /= spacing_[0]*spacing_[1];

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
			ival += -(val[0] - gridpos[1][0])*(val[1] - gridpos[1][1]) * gridval[2][dim];
			ival += -(val[0] - gridpos[2][0])*(val[1] - gridpos[2][1]) * gridval[1][dim];
			ival +=  (val[0] - gridpos[3][0])*(val[1] - gridpos[3][1]) * gridval[0][dim];

			ival /= spacing_[0]*spacing_[1];

			return ival;
		}

		~Grid2D(){}
	};	
}
