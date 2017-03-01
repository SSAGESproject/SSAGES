/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Julian Helfferich <julian.helfferich@gmail.com>
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

#include <exception>
#include <vector>

#include "Drivers/DriverException.h"
#include "JSON/Serializable.h"
#include "schema.h"
#include "Validator/ObjectRequirement.h"

// Forward declare.
namespace Json {
    class Value;
}

namespace SSAGES
{

//! Generic Grid.
/*!
 * \tparam T type of data to be stored in the Grid.
 *
 * A Grid is a general method to store data in SSAGES. It is used to discretize
 * a continuous number, typically a collective variable, into \c number_points
 * grid points. For each grid point an arbitrary type of data can be stored,
 * specified via the template parameter \c T.
 *
 * The Grid can be of arbitrary dimension. For each dimension, the lower bound,
 * the upper bound and the number of grid points need to be specified.
 * Furthermore, the grid can be defined as periodic or non-periodic in the
 * respective dimension. By default, the grid is non-periodic. The grid points
 * are indexed from 0 to number_points-1 following the standard C/C++
 * convention. The indices -1 and and \c number_points can be used to access
 * the underflow and overflow intervals (see below).
 *
 * The grid spacing \c Delta is given by (upper - lower)/number_points. Thus,
 * grid point \c n corresponds to the interval
 * [lower + n*Delta, lower + (n+1)*Delta) and the position of the grid point is
 * at the center of this interval, at lower + (n+0.5)*Delta. Note that n
 * follows the C/C++ convention, i.e. n = 0 for the first interval. The grid
 * indices pertaining to a given point can be obtained via Grid::GetIndices().
 *
 * In non-periodic dimensions, an overflow and an underflow interval exist. The
 * underflow interval corresponds to the interval (-infinity, lower), i.e. all
 * points below \c lower. Similarly, the overflow interval corresponds to the
 * interval [upper, infinity). The underflow grid point can be accessed via the
 * index -1, the overflow grid point via the index \c number_points.
 *
 * \ingroup Core
 */
template<typename T>
class Grid : public Serializable
{
private:
    //! Dimension of the grid
    const size_t dimension_;

    //! Number of points in each dimension.
    std::vector<int> numPoints_;

    //! Edges of the Grid in each dimension.
    std::pair< std::vector<double>,std::vector<double> > edges_;

    //! Periodicity of the Grid.
    std::vector<bool> isPeriodic_;

    //! Internal storage of the data
    std::vector<T> data_;

    //! Wrap the index around periodic boundaries
    std::vector<int> wrapIndices(const std::vector<int> &indices) const
    {
        std::vector<int> newIndices(dimension_);
        for (size_t i=0; i<dimension_; ++i) {
            if (!GetPeriodic(i)) {
                continue;
            }

            int index = indices.at(i);
            while (index < 0) { index += GetNumPoints(i); }
            while (index >= GetNumPoints(i)) { index -= GetNumPoints(i); }
            newIndices.at(i) = index;
        }

        return newIndices;
    }

    //! Map d-dimensional indices to 1-d data vector
    /*!
     * Map a set of indices to the index of the 1d data vector. Keep in mind,
     * that the data includes underflow (index -1) and overflow (index
     * numPoints) bins in periodic dimension.
     *
     * This function does not check if the indices are in bounds. This is done
     * in the function(s) calling mapTo1d().
     */
    size_t mapTo1d(const std::vector<int> &indices) const
    {
        size_t idx = 0;
        size_t fac = 1;
        for (size_t i=0; i<dimension_; ++i) {
            if (GetPeriodic(i)) {
                idx += indices.at(i) * fac;
                fac *= GetNumPoints(i);
            } else {
                idx += (indices.at(i) + 1) * fac;
                fac *= GetNumPoints(i) + 2;
            }
        }
        return idx;
    }

public:
    //! Constructor
    /*!
     * \param numPoints Number of grid points in each dimension.
     * \param lower Lower edges of the grid.
     * \param upper Upper edges of the grid.
     * \param isPeriodic Bools specifying the periodicity in the respective
     *                   dimension.
     *
     * The dimension of the grid is determined by the size of the parameter
     * vectors.
     */
    Grid(std::vector<int> numPoints,
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
        size_t data_size = 1;
        for (size_t d = 0; d < GetDimension(); ++d) {
            size_t storage_size = GetNumPoints(d);
            if (!GetPeriodic(d)) { storage_size += 2; }
            data_size *= storage_size;
        }

        data_.resize(data_size);
    }

    //! Get total number of grid points (including under- and overflow bins)
    /*!
     * This function returns the total number of grid points, which is identical
     * to the size of the internal data vector.
     */
    size_t size() const
    {
        return data_.size();
    }

    //! Get pointer to the internal data storage vector
    /*!
     * It is discouraged to directly access the internal data storage. It might,
     * however be necessary. For example when communicating the data over MPI.
     */
    T *data()
    {
        return data_.data();
    }

    //! Get pointer to const of the internal data storage vector
    /*!
     * It is discouraged to directly access the internal data storage. It might,
     * however be necessary. For example when communicating data over MPI.
     */
    T const *data() const
    {
        return data_.data();
    }

    //! Get the dimension.
    size_t GetDimension() const
    {
        return dimension_;
    }

    //! Get the number of points for all dimensions.
    const std::vector<int>& GetNumPoints() const
    {
        return numPoints_;
    }

    //! Get the number of points for a specific dimension.
    /*!
     * \param dim Index of the dimension.
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

        return numPoints_.at(dim);
    }

    //! Return the lower edges of the Grid.
    const std::vector<double>& GetLower() const
    {
        return edges_.first;
    }

    //! Get the lower edge for a specific dimension.
    /*!
     * \param dim Index of the dimension.
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
        return GetLower().at(dim);
    }

    //! Return the upper edges of the Grid.
    const std::vector<double>& GetUpper() const
    {
        return edges_.second;
    }


    //! Get the upper edge for a specific dimension.
    /*!
     * \param dim Index of the dimension.
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
        return GetUpper().at(dim);
    }

    //! Return the periodicity of the Grid.
    const std::vector<bool>& GetPeriodic() const
    {
        return isPeriodic_;
    }

    //! Get the periodicity in a specific dimension.
    /*!
     * \param dim Index of the dimension.
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
                if (xpos < GetLower(i)) {
                    indices.at(i) = -1;
                    continue;
                } else if (xpos > GetUpper(i)) {
                    indices.at(i) = GetNumPoints(i);
                    continue;
                }
            }

            // To make sure, the value is rounded in the correct direction.
            double round = 0.5;
            if (xpos < 0) { round = -0.5; }

            double spacing = (GetUpper(i) - GetLower(i)) / GetNumPoints(i);

            indices.at(i) = (xpos - GetLower(i)) / spacing + round;
        }

        return wrapIndices(indices);
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
        if (indices.size() != dimension_) {
            throw std::invalid_argument("Dimension of indices does not match "
                    "dimension of the grid.");
        }

        // Check if an index is out of bounds
        for (size_t i=0; i<dimension_; ++i) {
            int index = indices.at(i);
            if ( (isPeriodic_.at(i) && (index < 0 || index >= numPoints_.at(i))) ||
                 (!isPeriodic_.at(i) && (index < -1 || index > numPoints_.at(i))) )
            {
                throw std::out_of_range("Grid index out of range.");
            }
        }
        return data_.at(mapTo1d(indices));
    }

    //! Access Grid element read/write
    /*!
     * \param indices Vector of integers specifying the grid point.
     */
    T& at(const std::vector<int> &indices)
    {
        return const_cast<T&>(static_cast<const Grid<T>* >(this)->at(indices));
    }

    //! Const access of Grid element via initializer list
    /*!
     * \tparam R Datatype in the initializer list
     * \param x initializer list
     *
     * This function avoids abiguity if at() is called with a brace-enclosed
     * initializer list.
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
     *
     * This function avoids abiguity if at() is called with a brace-enclosed
     * initializer list.
     */
    template<typename R>
    T& at(std::initializer_list<R>&& x)
    {
        return at(static_cast<std::vector<R> >(x));
    }

    //! Access 1d Grid by index, read-only
    /*!
     * \param index Index specifying the grid point.
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
     */
    T& at(int index)
    {
        return const_cast<T&>(static_cast<const Grid<T>* >(this)->at(index));
    }

    //! Access Grid element pertaining to a specific point -- read-only
    /*!
     * \param x Vector of doubles specifying a point.
     *
     * This function is provided for convenience. It is identical to
     * Grid::at(Grid::GetIndices(x)).
     */
    const T& at(const std::vector<double> &x) const
    {
        return data_.at(GetIndices(x));
    }

    //! Access Grid element pertaining to a specific point -- read/write
    /*!
     * \param x Vector of doubles specifying a point.
     *
     * This function is provided for convenience. It is identical to
     * Grid::at(Grid::GetIndices(x)).
     */
    T& at(const std::vector<double> &x)
    {
        return data_.at(GetIndices(x));
    }

    //! Access 1d-Grid by point - read-only
    /*!
     * \param x Access grid point pertaining to this value.
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
     */
    T& at(double x)
    {
        return const_cast<T&>(static_cast<const Grid<T>* >(this)->at(x));
    }

    //! Access Grid element per [] read-only
    /*!
     * \param indices Vector of integers specifying the grid point.
     *
     * Example: grid[{0,1}]
     */
    const T& operator[](const std::vector<int> &indices) const
    {
        return at(indices);
    }

    //! Access Grid element per [] read-write
    /*!
     * \param indices Vector of integers specifying the grid point.
     *
     * Example: grid[{0,1}]
     */
    T& operator[](const std::vector<int> &indices)
    {
        return at(indices);
    }

    //! Const access of Grid element via initializer list
    /*!
     * \tparam R Datatype in the initializer list
     * \param x initializer list
     *
     * This function avoids abiguity if operator[] is called with a brace-enclosed
     * initializer list.
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
     *
     * This function avoids abiguity if operator[] is called with a brace-enclosed
     * initializer list.
     */
    template<typename R>
    T& operator[](std::initializer_list<R>&& x)
    {
        return at(static_cast<std::vector<R> >(x));
    }

    //! Access 1d-Grid per [] operator, read-only
    /*!
     * \param index Index of the grid point.
     */
    const T& operator[](int index) const
    {
        return at(index);
    }

    //! Access 1d-Grid per [] operator, read-write
    /*!
     * \param index Index of the grid point.
     */
    T& operator[](int index)
    {
        return at(index);
    }

    //! Access Grid element pertaining to a specific point per [] read-only
    /*!
     * \param indices Vector of integers specifying the grid point.
     *
     * Example: grid[{0.2,-1.43}]. Note, that you must not pass integers, as
     * integers will be interpreted as grid indices. Thus, call grid[{1.0, 2.0}]
     * and not grid[{1,2}].
     */
    const T& operator[](const std::vector<double> &x) const
    {
        return at(x);
    }

    //! Access Grid element pertaining to a specific point per [] read-write
    /*!
     * \param indices Vector of integers specifying the grid point.
     *
     * Example: grid[{0.2,-1.43}]. Note, that you must not pass integers, as
     * integers will be interpreted as grid indices. Thus, call grid[{1.0, 2.0}]
     * and not grid[{1,2}].
     */
    T& operator[](const std::vector<double> &x)
    {
        return at(x);
    }

    //! Access 1d-Grid via specific point, read-only
    /*!
     * \param x Point specifying the desired Grid point.
     */
    const T& operator[](double x) const
    {
        return at(x);
    }

    //! Access 1d-Grid via specific point, read-write
    /*!
     * \param x Point specifying the desired Grid point.
     */
    T& operator[](double x)
    {
        return at(x);
    }

    //! Set up the grid
    /*!
     * \param json JSON value containing all input information.
     *
     * This function builds a grid from a JSON node. It will return a nullptr
     * if an unknown error occured, but generally, it will throw a
     * BuildException of failure.
     */
    static Grid<T>* BuildGrid(const Json::Value& json)
    {
        return BuildGrid(json, "#/Grid");
    }

    //! Set up the grid
    /*!
     * \param json JSON Value containing all input information.
     * \param path Path for JSON path specification.
     *
     * This function builds a grid from a JSON node. It will return a nullptr
     * if an unknown error occured, but generally, it will throw a
     * BuildException on failure.
     */
    static Grid<T>* BuildGrid(const Json::Value& json, const std::string& path)
    {
        Json::ObjectRequirement validator;
        Json::Value schema;
        Json::Reader reader;

        reader.parse(JsonSchema::grid, schema);
        validator.Parse(schema, path);

        // Validate inputs.
        validator.Validate(json, path);
        if (validator.HasErrors()) {
            throw BuildException(validator.GetErrors());
        }

        // Read in Lower Grid edges.
        std::vector<double> lower;
        for (auto &lv : json["lower"]) {
           lower.push_back(lv.asDouble());
        }

        size_t dimension = lower.size();

        // Read in upper grid edges.
        std::vector<double> upper;
        for (auto &uv : json["upper"]) {
            upper.push_back(uv.asDouble());
        }

        if (upper.size() != dimension) {
            throw BuildException({"Number of upper values does not match "
                                  "number of lower values!"});
        }

        // Read in number of points.
        std::vector<int> number_points;
        for (auto &np : json["number_points"]) {
            number_points.push_back(np.asInt());
        }

        if (number_points.size() != dimension) {
            throw BuildException({"Arrays \"lower\" and \"number_points\" do "
                                  "not have the same size!"});
        }

        // Read in periodicity.
        std::vector<bool> isPeriodic;
        for (auto &periodic : json["periodic"]) {
            isPeriodic.push_back(periodic.asBool());
        }

        if (isPeriodic.size() == 0) {
            isPeriodic = std::vector<bool>(dimension, false);
        } else if (isPeriodic.size() != dimension) {
            throw BuildException({"Arrays \"lower\" and \"periodic\" do not "
                                  "have the same size!"});
        }

        // Construct the grid.
        Grid<T>* grid = new Grid(number_points, lower, upper, isPeriodic);

        return grid;
    }

    //! \copydoc Serializable::Serialize()
    /*!
     * \warning Serialization not yet implemented.
     */
    void Serialize(Json::Value& /*json*/) const override
    {

    }

    //! Return iterator at first element of internal storage
    typename std::vector<T>::iterator begin()
    {
        return data_.begin();
    }

    //! Return iterator after last element of internal storage
    typename std::vector<T>::iterator end()
    {
        return data_.end();
    }

    //! Return const iterator at first element of internal storage
    typename std::vector<T>::const_iterator begin() const
    {
        return data_.begin();
    }

    //! Return const iterator after last element of internal storage
    typename std::vector<T>::const_iterator end() const
    {
        return data_.end();
    }
};

} // End namespace SSAGES
