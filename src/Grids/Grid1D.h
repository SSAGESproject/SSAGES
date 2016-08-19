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
			_NDim = 1;
			for(size_t i =0; i <lower.size(); i++)
			{
				_lower.push_back(lower[i]);
				_upper.push_back(upper[i]);
				_num_points.push_back(num_points[i]);
				_spacing.push_back((_upper[i] - _lower[i])/double(_num_points[i] - 1));
				_periodic.push_back(false);
			}

			// Construct flat vector 
			_flatvector.resize(_num_points[0]);

			for(int i = 0; i < _num_points[0]; i++)
			{
				_flatvector[i].first = std::vector<double>(2, 0.0);
				_flatvector[i].second.push_back(_lower[0] + _spacing[0]*i);
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
