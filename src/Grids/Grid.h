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

#include "../JSON/Serializable.h"
#include <vector>
#include <iostream>
#include <stdexcept>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	//! Calculate array index from n-Dimensional Grid indices.
	/*!
	 * \param indices Vector of indices.
	 * \param num_points Number of grid points in each dimension.
	 *
	 * \return Index of 1d array storing data.
	 */
	inline int FlattenIndices(std::vector<int> indices, std::vector<int> num_points)
	{

		int loci = 0;
		for(size_t i = 0; i < indices.size(); i++)
		{
			if(indices[i]>=num_points[i])
				throw std::out_of_range("Grid indices out of range");
			
			int locj = indices[i];
			for(size_t j = i+1;j<indices.size();j++)
				locj *= num_points[j];
			loci += locj;
		}
		
		return loci;
	}

	//! Generic Grid.
	/*!
	 * \ingroup Core
	 */
	class Grid: public Serializable
	{

	protected:

		std::vector<double> lower_; //!< Lower edge of the grid.
		std::vector<double> upper_; //!< Upper edge of the grid.
		std::vector<bool> periodic_; //!< Is the grid periodic in the corresponding dimension?
		std::vector<int> num_points_; //!< Number of grid points.
		std::vector<double> spacing_; //!< Grid spacing.
		int NDim_; //!< Grid dimension.

		//! Array storing grid data.
		std::vector<std::pair<std::vector<double>, std::vector<double>>> flatvector_;
	
	public:
		//! Iterator for traversing the Grid.
		using const_iterator = std::vector<std::pair<std::vector<double>, std::vector<double>>>::const_iterator;

		//! Destructor.
		virtual ~Grid(){}

		//! Return the nearest indices for given values.
		/*!
		 * \param val N-dimensional value.
		 * \return N-dimensional vector with nearest grid indices.
		 *
		 * Returns the grid point which is the nearest to a given value. The
		 * value must be given as a N-dimensional vector and the grid point is
		 * returned as an N-dimensional vector. Here, N is the dimension of the
		 * Grid.
		 */
		std::vector<int> GetIndices(const std::vector<double> &val)
		{
			std::vector<int> vertices;

			for(size_t i = 0; i < val.size(); i++)
			{
				int vertex = 0;
				if(val[i]<lower_[i])
					vertex = int((val[i] - lower_[i])/spacing_[i] - 0.5);
				else
					vertex = int((val[i] - lower_[i])/spacing_[i] + 0.5);

				if(periodic_[i])
				{
					vertex = vertex % num_points_[i];
					if(vertex<0)
						vertex += num_points_[i];
				}

				if(vertex < 0) // out of bounds
					vertex = 0;
				else if(vertex >= num_points_[i]) // out of bounds
					vertex = num_points_[i] - 1;

				vertices.push_back(vertex);
			}

			return vertices;
		}

		//! Get the location at the current indices.
		/*!
		 * \param indices Indices specifying grid point.
		 * \return vector of positions at the specified grid point.
		 */
		std::vector<double> GetLocation(const std::vector<int> &indices) const
		{
			return flatvector_[FlattenIndices(indices,num_points_)].second;
		}

		//! Get the value at the current indices.
		/*!
		 * \param indices Indices specifying grid point.
		 * \return Value at the specified grid point.
		 */
		double GetValue(const std::vector<int>& indices) const
		{
			return flatvector_[FlattenIndices(indices,num_points_)].first[0];
		}

		//! Set the value at the current incices.
		/*!
		 * \param indices Indices specifying the grid point.
		 * \param value New value for the grid point.
		 */
		void SetValue(const std::vector<int>& indices, double value)
		{
			flatvector_[FlattenIndices(indices,num_points_)].first[0] = value;
		}

		//! Get the extra at the current indices.
		/*!
		 * \param indices Indices specifying grid point.
		 * \return Value at the specified grid point.
		 */
		std::vector<double> GetExtra(const std::vector<int>& indices) const
		{
			const std::vector<double>& temp = flatvector_[FlattenIndices(indices,num_points_)].first;
			std::vector<double> temp2(&temp[1], &temp[temp.size()]);
			return temp2;
		}

		//! Set the extra at the current incices.
		/*!
		 * \param indices Indices specifying the grid point.
		 * \param value New value for the grid point.
		 */
		void SetExtra(const std::vector<int>& indices, std::vector<double> value)
		{
			int flatten = FlattenIndices(indices, num_points_);
			if(value.size() != flatvector_[flatten].first.size()-1)
				throw std::out_of_range("Vector field in grid is not the same size as value");
			for(size_t i = 0; i < value.size(); i++)
				flatvector_[flatten].first[i+1] = value[i];
		}

		//! Set the entire grid given two vectors of values and dimension.
		/*!
		 * \param first_values Values specifying the grid point value and derivatives/vector field.
		 */
		void SetGrid(const std::vector<double>& first_values)
		{
			std::cout<<first_values.size()/(NDim_+1)<<" "<<flatvector_.size()<<std::endl;
			if(first_values.size()/(NDim_+1) != flatvector_.size())
				throw std::out_of_range("Attempting to set grid with more/less values than grid has.");

			for(size_t i = 0; i < flatvector_.size(); i++)
				for(int j = 0; j < NDim_ + 1; j++)
					flatvector_[i].first[j] = first_values[i*(NDim_+1) + j];
		}

		//! Set the entire grid given two vectors of values and dimension.
		/*!
		 * \param periodic_values Values specifying the grid PBC.
		 */
		void SetPeriodic(const std::vector<bool>& periodic_values)
		{
			if(periodic_values.size() != periodic_.size())
				throw std::out_of_range("Value dimension used to set periodic boundaries are not equal.");

			for(size_t i = 0; i < periodic_values.size(); i++)
				periodic_[i] = periodic_values[i];
		}

		//! Write Grid to console
		/*!
		 * Currently only for debugging
		 */
		void PrintGrid() const
		{
			for(size_t i = 0; i<flatvector_.size(); i++)
			{
				for(size_t j = 0; j<flatvector_[i].first.size(); j++)
					std::cout << flatvector_[i].first[j]<<" ";
				for(size_t j = 0; j<flatvector_[i].second.size(); j++)
					std::cout<<flatvector_[i].second[j]<< " ";
				std::cout<<std::endl;
			}
		}

		//! Return lower edges of the Grid.
		/*!
		 * \return Vector containing the lower edges of the grid in each dimension.
		 */
		std::vector<double> GetLower()
		{
			return lower_;
		}

		//! Return upper edges of the Grid.
		/*!
		 * \return Vector containing the upper edges of the grid in each dimension.
		 */
		std::vector<double> GetUpper()
		{
			return upper_;
		}

		//! Get periodicity of the Grid.
		/*!
		 * \return Vector with boolean values.
		 *
		 * The vector contains \c True if the Grid is periodic in the
		 * corresponding dimension.
		 */
		std::vector<bool> GetPeriodic()
		{
			return periodic_;
		}

		//! Get number of Grid points.
		/*!
		 * \return Vector containing number of grid points in each dimension.
		 */
		std::vector<int> GetNumPoints()
		{
			return num_points_;
		}

		//! Get grid spacing.
		/*!
		 * \return Vector containing grid spacing in each dimension.
		 */
		std::vector<double> GetSpacing()
		{
			return spacing_;
		}

		//! Get dimensionality of the Grid.
		/*!
		 * \return Grid dimension.
		 */
		int GetDimension()
		{
			return NDim_;
		}

		//! Set up the Grid.
		/*!
		 * \param json JSON input value.
		 * \return Pointer to the Grid. \c nullptr if unknown error occured.
		 *
		 * This function builds a Grid from a JSON node. If an unknown error
		 * occured, the return value is \c nullptr, but generally, the function
		 * will throw a BuildException on failure.
		 * \note Object lifetime is the caller's responsibility.
		 */
		static Grid* BuildGrid(const Json::Value& json);

		//! Set up Grid.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification
		 * \return Pointer the the Grid built.
		 *
		 * This function overloads Grid::BuildGrid(const Json::Value&) to allow
		 * JSON path specification.
		 */
		static Grid* BuildGrid(const Json::Value& json, 
							   const std::string& path);

		//! Get a voxel from the Grid.
		/*!
		 * \param val Position value.
		 * \return Voxel
		 *
		 * This function returns a voxel from the Grid. For each dimension, the
		 * voxel contains the next smaller and the next larger index based on the
		 * given position. In other words, the position is inside the
		 * (hyper-)cube defined by the voxel grid points.
		 */
		virtual std::vector<std::vector<int>> GetVoxel(const std::vector<double> &val) const = 0;

		//! Interpolate value.
		/*!
		 * \param val Position.
		 *
		 * \return Interpolated value.
		 *
		 * This function returns the interpolated value at the given position.
		 * The value is a linear interpolation between the value at the next
		 * smaller and the next larger grid index, taking all dimensions of the
		 * grid into account.
		 */
		virtual double InterpolateValue(const std::vector<double> &val) const = 0;

		//! Interpolate derivative value.
		/*!
		 * \param val Position.
		 * \param dim Dimension in which respect to calculate the derivative.
		 *
		 * \return Derivative of the interpolated value with respect to the
		 *         given dimension.
		 *
		 * This function determines the derivative of the interpolated value with
		 * respect to a given dimension.
		 *
		 * \sa InterpolateValue()
		 */
		virtual double InterpolateDeriv(const std::vector<double> &val, int dim) const = 0;

		//! \copydoc Serializable::Serialize()
		virtual void Serialize(Json::Value& json) const override
		{
			for(auto& p : lower_)
				json["lower"].append(p);

			for(auto& p : upper_)
				json["upper"].append(p);

			for(auto& p : num_points_)
				json["number_points"].append(p);

			for(auto p : periodic_)
				json["periodic"].append(p);

			for(auto& point : flatvector_)
				for(auto& p : point.first)
					json["values"].append(p);

		}

		//! Get iterator to first element of the Grid.
		/*!
		 * \return Iterator to first element in the Grid.
		 */
		const_iterator begin() const { return flatvector_.begin();}

		//! Get iterator one beyond last element of the Grid.
		/*!
		 * \return Iterator pointing to after the last element in the Grid.
		 */
		const_iterator end() const { return flatvector_.end(); }
	};
	
}
