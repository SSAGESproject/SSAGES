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
		 * \param num_points List of numbers of Grid points.
		 */
		Grid1D(std::vector<double> lower, std::vector<double> upper, std::vector<int> num_points)
		{
			NDim_ = 1;
			for(size_t i =0; i <lower.size(); i++)
			{
				lower_.push_back(lower[i]);
				upper_.push_back(upper[i]);
				num_points_.push_back(num_points[i]);
				spacing_.push_back((upper_[i] - lower_[i])/double(num_points_[i] - 1));
				periodic_.push_back(false);
			}

			// Construct flat vector 
			flatvector_.resize(num_points_[0]);

			for(int i = 0; i < num_points_[0]; i++)
			{
				flatvector_[i].first = std::vector<double>(2, 0.0);
				flatvector_[i].second.push_back(lower_[0] + spacing_[0]*i);
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
			voxel.push_back(vertices);

			for(size_t v = 0; v < vertices.size(); v++)
			{
				vertices[v]++;
				voxel.push_back(vertices);
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

			//now, do 1d interpolation
			double ival = ((val[0]-gridpos[0][0])*gridval[1] +
			  (gridpos[1][0]-val[0])*gridval[0]) /
			      (gridpos[1][0]-gridpos[0][0]);

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

			//now, do 1d interpolation
			double ival = ((val[0]-gridpos[0][0])*gridval[1][dim] +
			  (gridpos[1][0]-val[0])*gridval[0][dim]) /
			      (gridpos[1][0]-gridpos[0][0]);

			return ival;
		}

		~Grid1D(){}
	};	
}
