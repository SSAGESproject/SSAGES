/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Julian Helfferich <julian.helfferich@gmail.com>
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

#include <cmath>
#include <exception>
#include <vector>

namespace SSAGES
{

//! Base class for Grids
/*!
 * \tparam Type of data stored on the grid
 *
 * Base class for all grids. Currently, these are 'Grid' and 'Histogram'.
 *
 * \ingroup Core
 */
template<typename T>
class GridBase
{
protected:
    //! Internal storage of the data
    std::vector<T> data_;

    //! Dimension of the grid
    size_t dimension_;

    //! Number of points in each dimension.
    std::vector<int> numPoints_;

    //! Edges of the Grid in each dimension.
    std::pair< std::vector<double>,std::vector<double> > edges_;

    //! Periodicity of the Grid.
    std::vector<bool> isPeriodic_;

    //! Wrap the index around periodic boundaries
    std::vector<int> wrapIndices(const std::vector<int> &indices) const
    {
        std::vector<int> newIndices(indices);
        for (size_t i=0; i<dimension_; ++i) {
            if (!GetPeriodic(i)) {
                continue;
            }

            int index = indices.at(i);
            while (index < 0) { index += numPoints_[i]; }
            while (index >= numPoints_[i]) { index -= numPoints_[i]; }
			if(index == numPoints_[i]-1)
				index = 0;
            newIndices.at(i) = index;
        }

        return newIndices;
    }

    //! This function needs to be implemented by child classes.
    /*!
     * \param indices The indices specifying the grid point.
     * \return Index of the grid point in the 1d storage vector.
     *
     * mapTo1d maps the indices onto a single index to access the data stored
     * in the data_ vector. This mapping will be different for different grid
     * implementations.
     */
    virtual size_t mapTo1d(const std::vector<int> &indices) const = 0;

protected:
    //! Constructor
    /*!
     * \param numPoints Number of grid points in each dimension.
     * \param lower Lower edges of the grid.
     * \param upper Upper edges of the grid.
     * \param isPeriodic Bools specifying the periodicity in the respective
     *                   dimension.
     *
     * The constructor is protected by design. This makes sure that only child
     * classes of GridBase are constructed.
     */
    GridBase(std::vector<int> numPoints,
             std::vector<double> lower,
             std::vector<double> upper,
             std::vector<bool> isPeriodic)
        : dimension_(numPoints.size()),
          numPoints_(numPoints),
          edges_(std::pair< std::vector<double>, std::vector<double> >(lower, upper)),
          isPeriodic_(isPeriodic)
    {
        // Check that vector sizes are correct
        if (edges_.first.size() != dimension_ ||
            edges_.second.size() != dimension_) {
            throw std::invalid_argument("Size of vector containing upper or "
                "lower edges, does not match size of vector containing "
                "number of grid points.");
        }
        if (isPeriodic_.size() == 0) {
            // Default: Non-periodic in all dimensions
            isPeriodic.resize(dimension_, false);
        } else if (isPeriodic_.size() != dimension_) {
            throw std::invalid_argument("Size of vector isPeriodic does not "
                    "match size of vector containing number of grid points.");
        }
		for(size_t i=0 ; i < isPeriodic_.size() ; i++)
		{
			if(isPeriodic_[i])
			{
				double spacing = (edges_.second[i] - edges_.first[i]) / numPoints_[i];
				numPoints_[i]++;
				edges_.first[i] -= spacing/2;
				edges_.second[i] += spacing/2;
			}
		}
		
    }

public:

	
	void syncGrid()
	{
		double dim = this->GetDimension();
		
		if(dim == 1 && isPeriodic_[0])
		{
			this->at(numPoints_[0]-1) = this->at(0);
		}
		
		else if(dim == 1)
		{
			return;			
		}
		
		double loopsize = 1;
		std::vector<int> navigate(dim);		
		for(size_t i = 0; i<dim ; ++i)
		{
			if(!isPeriodic_[i])
				continue;
				
			loopsize = 1;
			for(size_t j = 0; j<dim; ++j)
			{
				if(j!=i)
					loopsize *= numPoints_[j];
			}
			std::fill(navigate.begin(), navigate.end(), 0);
			
			double div = 1;			
			for(size_t j = 0; j<loopsize; ++j)
			{
				div = 1;
				if(i==0)
				{
					navigate[1] = j%(numPoints_[1]);
				}
				else
				{
					navigate[0] = j%(numPoints_[0]);
				}
				for(size_t l = 1; l<dim; ++l)
				{
					if((i==0)&&(l==1))
					{
						continue;
					}
					if( l!=i )
					{
						div *= numPoints_[l];	
						navigate[l] = loopsize/div;
					}
				}
				navigate[i] = 0;
				auto temp = this->at(navigate);
				navigate[i] = numPoints_[i]-1;
				this->at(navigate) = temp;
			}
		}			
	}	
	
	 
	//! Get the dimension.
    /*!
     * \return Dimensionality of the grid.
     */
    size_t GetDimension() const
    {
        return dimension_;
    }

    //! Get the number of points for all dimensions.
    /*!
     * \return Vector of ints containing the number of grid points for each
     *         dimension.
     */
    const std::vector<int> GetNumPoints() const
    {
		std::vector<int> numPoints = numPoints_;
		for(size_t i = 0; i<dimension_;++i)
			if(GetPeriodic(i)){numPoints[i]--;}
        return numPoints;
    }

    //! Get the number of points for a specific dimension.
    /*!
     * \param dim Index of the dimension.
     * \return Number of grid points in the requested dimension.
     *
     * \note The first dimension uses the index 0.
     */
    int GetNumPoints(size_t dim) const
    {
        if (dim >= GetDimension()) {
            std::cerr << "Warning! Grid size requested for a dimension larger "
                         "than the grid dimensionality!\n";
            return 0;
        }
		if(GetPeriodic(dim)){return numPoints_[dim]-1;}
		else
			return numPoints_.at(dim);
    }

    //! Return the lower edges of the Grid.
    /*!
     * \return Vector containing the lower edges of the grid.
     */
    const std::vector<double> GetLower() const
    {
		std::vector<double> lower(dimension_);
		for(size_t i = 0; i<dimension_;++i) 
			if(GetPeriodic(i)) {lower[i] = edges_.first[i] - ((edges_.first[i] - edges_.second[i]) / numPoints_[i])/2;}
			else {lower[i] = edges_.first[i];}
        return lower;
    }

    //! Get the lower edge for a specific dimension.
    /*!
     * \param dim Index of the dimension.
     * \return Value of the lower edge in the requested dimension.
     *
     * \note The first dimension has the index 0.
     */
    double GetLower(size_t dim) const
    {
        if (dim >= GetDimension()) {
            std::cerr << "Warning! Lower edge requested for a dimension larger "
                         "than the grid dimensionality!\n";
            return 0.0;
        }
        if(GetPeriodic(dim)){return edges_.first[dim] - ((edges_.first[dim] - edges_.second[dim]) / numPoints_[dim])/2;}
		else{return edges_.first[dim];}
    }

    //! Return the upper edges of the Grid.
    /*!
     * \return Vector containing the upper edges of the grid.
     */
    const std::vector<double> GetUpper() const
    {
		std::vector<double> upper(dimension_);
		for(size_t i = 0; i<dimension_;++i) 
			if(GetPeriodic(i)) {upper[i] = edges_.second[i] + ((edges_.first[i] - edges_.second[i]) / numPoints_[i])/2;}
			else {upper[i] = edges_.second[i];}
        return upper;
    }


    //! Get the upper edge for a specific dimension.
    /*!
     * \param dim Index of the dimension.
     * \return Value of the upper edge in the given dimension.
     *
     * \note The dimensions are indexed starting with 0.
     */
    double GetUpper(size_t dim) const
    {
        if (dim >= GetDimension()) {
            std::cerr << "Warning! Upper edge requested for a dimension larger "
                         "than the grid dimensionality!\n";
            return 0.0;
        }
		if(GetPeriodic(dim)){return edges_.second[dim] + ((edges_.first[dim] - edges_.second[dim]) / numPoints_[dim])/2;}
		else{return edges_.second[dim];}
    }

    //! Return the periodicity of the Grid.
    /*!
     * \return Vector of bools. The values are \c True (\c False ) if the grid
     *         is periodic (non-periodic) in the given dimension.
     */
    const std::vector<bool>& GetPeriodic() const
    {
        return isPeriodic_;
    }

    //! Get the periodicity in a specific dimension.
    /*!
     * \param dim Index of the dimension.
     * \return \c True (\c False ) if the grid is periodic (non-periodic) in
     *         the specified dimension.
     *
     * \note The dimensions are indexed starting with 0.
     */
    bool GetPeriodic(size_t dim) const
    {
        if (dim >= GetDimension()) {
            std::cerr << "Warning! Periodicity requested for a dimension larger "
                         "than the grid dimensionality!\n";
            return false;
        }
        return GetPeriodic().at(dim);
    }

    //! Get the size of the internal storage vector
    /*!
     * \return Size of the internal storage vector.
     *
     * This function returns the size of the internal storage vector. This is
     * also the total number of grid points including the over/underflow bins
     * in case of a histogram.
     */
    size_t size() const
    {
        return data_.size();
    }

    //! Get pointer to the internal data storage vector
    /*!
     * \return Pointer to data in the internal storage vector.
     *
     * It is discouraged to directly access the internal data storage. It might,
     * however be necessary. For example when communicating the data over MPI.
     */
    T *data()
    {
        return data_.data();
    }

    //! Get pointer to const of the internal data storage vector
    /*!
     * \return Const pointer to data in the internal storage vector.
     *
     * It is discouraged to directly access the internal data storage. It might,
     * however be necessary. For example when communicating data over MPI.
     */
    T const *data() const
    {
        return data_.data();
    }

    //! Return the Grid indices for a given point.
    /*!
     * \param x Point in space.
     * \return Indices of the grid point to which the point in space pertains.
     *
     * The grid discretizes the continuous space. For a given point in this
     * continuous space, this function will return the indices of the grid point
     * covering the point in space.
     *
     * If the grid is non-periodic in a given dimension and x is lower than the
     * lower edge in this dimension, the function will return -1, the index of
     * the underflow bin. Similarly, it will return numPoints, the index of the
     * overflow bin, if x is larger than the upper edge.
     *
     * In periodic dimensions, the index can take any integer value and will be
     * wrapped to the interval [0, numPoints).
     */
    std::vector<int> GetIndices(const std::vector<double> &x) const
    {
        // Check that input vector has the correct dimensionality
        if (x.size() != dimension_) {
            throw std::invalid_argument("Specified point has a larger "
                    "dimensionality than the grid.");
        }

        std::vector<int> indices(dimension_);
        for (size_t i = 0; i < dimension_; ++i) {
            double xpos = x.at(i);
            if (!GetPeriodic(i)) {
                if (xpos < edges_.first[i]) {
                    indices.at(i) = -1;
                    continue;
                } else if (xpos > edges_.second[i]) {
                    indices.at(i) = numPoints_[i];
                    continue;
                }
            }

            // To make sure, the value is rounded in the correct direction.
            double spacing = (edges_.second[i] - edges_.first[i]) / numPoints_[i];
            indices.at(i) = std::floor( (xpos - edges_.first[i]) / spacing);
        }
		
        return wrapIndices(indices);
    }
	
	//! Return linear interpolation on a coordinate.
	/*!
	 * \param x Point in space.
     * \return Linearly interpolated value of grid at that point.
     *  This function performs a n-linear interpolation on the grid.
	 *  The formula is bilinear/trilinear interpolation generalized
     *  to N-dimensional space. 
     */
	double GetInterpolated(const std::vector<double> &x)
    {
		double interpvalue = 0;
		int tempindex;
		for (size_t i = 0 ; i < pow(2,dimension_) ; ++i)
		{
			
			std::vector<double> shiftedvector = x;
			tempindex = i;
				
			double accumulatedweight = 1;
			for(size_t j = 0; j < dimension_ ; ++j)
			{
				double spacing = (edges_.second[j] - edges_.first[j]) / numPoints_[j];
				shiftedvector[j] += ((tempindex%2)-0.5)*spacing;
				tempindex = tempindex/2;
			}
			
			std::vector<int> shiftedindices = GetIndices(shiftedvector);
			std::vector<double> shiftedcenters = GetCoordinates(shiftedindices);
			
			for(size_t j = 0; j < dimension_ ; ++j)
			{	
				double spacing = (edges_.second[j] - edges_.first[j]) / numPoints_[j];
				//Handle Edges
				if(shiftedindices[j] == -1)
					{
					accumulatedweight *= ((std::abs(x[j]-GetCoordinates(GetIndices(x))[j])/spacing));
					shiftedvector[j] += spacing/2;
					}
				else if (shiftedindices[j] == numPoints_[j])
					{
					accumulatedweight *= ((std::abs(x[j]-GetCoordinates(GetIndices(x))[j])/spacing));
					shiftedvector[j] -= spacing/2;
					}
				else
					{
					// Handle Periodicity
					if(std::abs(x[j]-shiftedcenters[j]) > spacing)
						accumulatedweight *= (1-(std::abs(std::abs(x[j]-shiftedcenters[j]) - (edges_.second[j] - edges_.first[j]))/spacing));
					else
						accumulatedweight *= (1-(std::abs(x[j]-shiftedcenters[j])/spacing));
					}
			}			
			interpvalue += accumulatedweight*at(shiftedvector);
		}
		return interpvalue;
    }

    //! Return the Grid index for a one-dimensional grid.
    /*!
     * \param x Point in space.
     * \return Grid index to which the point pertains.
     *
     * Return the Grid index pertaining to the given point in space. This
     * function is for convenience when accessing 1d-Grids. For
     * higher-dimensional grids, x needs to be a vector of doubles.
     */
    int GetIndex(double x) const
    {
        if (dimension_ != 1) {
            throw std::invalid_argument("1d Grid index can only be accessed for "
                   "1d-Grids can be accessed with a.");
        }
        return GetIndices({x}).at(0);
    }

    //! Return coordinates of the grid center points
    /*!
     * \param indices Grid indices specifying a grid point.
     * \return Vector of double specifying the position of the grid point.
     *
     * The grid is a discretization of real or cv space. Thus, each grid point
     * is associated with an interval of the underlying space. This function
     * returns the center point of this interval.
     */
    std::vector<double> GetCoordinates(const std::vector<int> &indices)
    {
        if (indices.size() != dimension_) {
            throw std::invalid_argument(
                    "Grid indices specified for GetCoordinates() do not have "
                    "the same dimensionality as the grid.");
        }

        std::vector<double> v(dimension_);

        for (size_t i = 0; i < dimension_; ++i) {
            double spacing = (edges_.second[i] - edges_.first[i]) / numPoints_[i];
            v.at(i) = edges_.first[i] + (indices[i] + 0.5)*spacing;
        }

        return v;
    }

    //! Return center point of 1d-grid
    /*!
     * \param index Index of the 1d grid.
     * \return Coordinate in real/CV space.
     *
     * \note This function is only available for 1d grids.
     */
    double GetCoordinate(int index)
    {
        return GetCoordinates({index}).at(0);
    }

    //! Access Grid element read-only
    /*!
     * \param indices Vector of integers specifying the grid point.
     * \return const reference of the value stored at the given grid point.
     *
     * In non-periodic dimensions, the index needs to be in the interval
     * [-1, numPoints]. Grid::at(-1) accessed the underflow bin,
     * Grid::at(numPoints) accesses the overflow bin.
     *
     * In periodic dimensions, the index may take any integer value and will be
     * mapped back to the interval [0, numPoints-1]. Thus, Grid::at(-1) will
     * access the same value as Grid::at(numPoints-1).
     */
    const T& at(const std::vector<int> &indices) const
    {
        // Check that indices are in bound.
        if (indices.size() != GetDimension()) {
            throw std::invalid_argument("Dimension of indices does not match "
                    "dimension of the grid.");
        }

        return data_.at(mapTo1d(indices));
    }

    //! Access Grid element read/write
    /*!
     * \param indices Vector of integers specifying the grid point.
     * \return Reference to value at the specified grid point.
     */
    T& at(const std::vector<int> &indices)
    {
        return const_cast<T&>(static_cast<const GridBase<T>* >(this)->at(indices));
    }

    //! Const access of Grid element via initializer list
    /*!
     * \tparam R Datatype in the initializer list
     * \param x initializer list
     * \return Const reference to value at the specified point.
     *
     * This function avoids abiguity if at() is called with a brace-enclosed
     * initializer list. The template parameter makes sure that this function
     * can be called with either ints, specifying a grid point, or doubles,
     * specifying coordinates in space, inside the initializer list.
     */
    template<typename R>
    const T& at(std::initializer_list<R>&& x) const
    {
        return at(static_cast<std::vector<R> >(x));
    }

    //! Access Grid element via initializer list
    /*!
     * \tparam R Datatype in the initializer list
     * \param x initializer list
     * \return Reference to value at the specified point.
     *
     * This function avoids abiguity if at() is called with a brace-enclosed
     * initializer list. The template parameter makes sure that this function
     * can be called with either ints, specifying a grid point, or doubles,
     * specifying coordinates in space, inside the initializer list.
     */
    template<typename R>
    T& at(std::initializer_list<R>&& x)
    {
        return at(static_cast<std::vector<R> >(x));
    }

    //! Access 1d Grid by index, read-only
    /*!
     * \param index Index specifying the grid point.
     * \return Const reference of value at the given index.
     *
     * \note This function can only be used for 1d-Grids.
     */
    const T& at(int index) const
    {
        if (dimension_ != 1) {
            throw std::invalid_argument("Only 1d-Grids can be accessed with a "
                    "single integer as the index.");
        }
        return at({index});
    }

    //! Access 1d Grid by index, read-write
    /*!
     * \param index Index specifying the grid point.
     * \return Reference of value at the given grid point.
     *
     * \note This function can only be used for 1d-Grids.
     */
    T& at(int index)
    {
        return const_cast<T&>(static_cast<const GridBase<T>* >(this)->at(index));
    }

    //! Access Grid element pertaining to a specific point -- read-only
    /*!
     * \param x Vector of doubles specifying a point.
     * \return Const reference of the value at the given coordinates.
     *
     * This function is provided for convenience. It is identical to
     * GridBase::at(GridBase::GetIndices(x)).
     */
    const T& at(const std::vector<double> &x) const
    {
        return at(GetIndices(x));
    }

    //! Access Grid element pertaining to a specific point -- read/write
    /*!
     * \param x Vector of doubles specifying a point.
     * \return Reference to the value at the given coordinates.
     *
     * This function is provided for convenience. It is identical to
     * GridBase::at(GridBase::GetIndices(x)).
     */
    T& at(const std::vector<double> &x)
    {
        return at(GetIndices(x));
    }

    //! Access 1d-Grid by point - read-only
    /*!
     * \param x Access grid point pertaining to this value.
     * \return Const reference to the value pertaining to the specified
     *         coordinate.
     *
     * \note This function can only be used for 1d-Grids.
     */
    const T& at(double x) const
    {
        if (dimension_ != 1) {
            throw std::invalid_argument("Only 1d-Grids can be accessed with a "
                    "single float as the specified point.");
        }
        return at({x});
    }

    //! Access 1d-Grid by point - read-write
    /*!
     * \param x Access grid point pertaining to this value.
     * \return Reference to the value pertaining to the specified coordinate.
     *
     * \note This function can only be used for 1d Grids.
     */
    T& at(double x)
    {
        return const_cast<T&>(static_cast<const GridBase<T>* >(this)->at(x));
    }

    //! Access Grid element per [] read-only
    /*!
     * \param indices Vector of integers specifying the grid point.
     * \return Const reference to value to the given grid point.
     */
    const T& operator[](const std::vector<int> &indices) const
    {
        return at(indices);
    }

    //! Access Grid element per [] read-write
    /*!
     * \param indices Vector of integers specifying the grid point.
     * \return Reference to value at the specified grid point.
     */
    T& operator[](const std::vector<int> &indices)
    {
        return at(indices);
    }

    //! Const access of Grid element via initializer list
    /*!
     * \tparam R Datatype in the initializer list
     * \param x initializer list
     * \return Const reference to the value at the specified point.
     *
     * This function avoids abiguity if operator[] is called with a brace-enclosed
     * initializer list.
     *
     * Example: grid[{0,1}] or grid[{-1.23, 4.2, 0.0}]
     */
    template<typename R>
    const T& operator[](std::initializer_list<R>&& x) const
    {
        return at(static_cast<std::vector<R> >(x));
    }

    //! Access Grid element via initializer list
    /*!
     * \tparam R Datatype in the initializer list
     * \param x initializer list
     * \return Reference to the value at the specified point.
     *
     * This function avoids abiguity if operator[] is called with a brace-enclosed
     * initializer list.
     *
     * Example: grid[{0,1}]
     */
    template<typename R>
    T& operator[](std::initializer_list<R>&& x)
    {
        return at(static_cast<std::vector<R> >(x));
    }

    //! Access 1d-Grid per [] operator, read-only
    /*!
     * \param index Index of the grid point.
     * \return Const reference to the value at the given grid point.
     *
     * \note This function can only be used for 1d grids.
     */
    const T& operator[](int index) const
    {
        return at(index);
    }

    //! Access 1d-Grid per [] operator, read-write
    /*!
     * \param index Index of the grid point.
     * \return Reference to the value at the given grid point.
     *
     * \note This function can only be used for 1d grids.
     */
    T& operator[](int index)
    {
        return at(index);
    }

    //! Access Grid element pertaining to a specific point per [] read-only
    /*!
     * \param x Vector of doubles specifying the point in space.
     * \return Const reference to the value pertaining to the given coordinates.
     */
    const T& operator[](const std::vector<double> &x) const
    {
        return at(x);
    }

    //! Access Grid element pertaining to a specific point per [] read-write
    /*!
     * \param x Vector of doubles specifying the point in space.
     * \return Reference to the value pertaining to the given coordinates.
     */
    T& operator[](const std::vector<double> &x)
    {
        return at(x);
    }

    //! Access 1d-Grid via specific point, read-only
    /*!
     * \param x Point specifying the desired Grid point.
     * \return Const reference to the value at the specified coordinate.
     *
     * \note This function can only be used for 1d grids.
     */
    const T& operator[](double x) const
    {
        return at(x);
    }

    //! Access 1d-Grid via specific point, read-write
    /*!
     * \param x Point specifying the desired Grid point.
     * \return Reference to value at the specified coordinate.
     *
     * \note This function can only be used for 1d grids.
     */
    T& operator[](double x)
    {
        return at(x);
    }
};

} // End namespace SSAGES
