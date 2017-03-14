/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016-2017 Julian Helfferich <julian.helfferich@gmail.com>
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

#include "Drivers/DriverException.h"
#include "GridBase.h"
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
class Grid : public GridBase<T>
{
private:
    //! Map d-dimensional indices to 1-d data vector
    /*!
     * Map a set of indices to the index of the 1d data vector. Keep in mind,
     * that the data includes underflow (index -1) and overflow (index
     * numPoints) bins in periodic dimension.
     *
     * This function does not check if the indices are in bounds. This is done
     * in the function(s) calling mapTo1d().
     */
    size_t mapTo1d(const std::vector<int> &indices) const override
    {
        // Check if an index is out of bounds
        for (size_t i=0; i < GridBase<T>::GetDimension(); ++i) {
            int index = indices.at(i);
            bool periodic = GridBase<T>::GetPeriodic(i);
            int numpoints = GridBase<T>::GetNumPoints(i);
            if ( (periodic && (index < 0 || index >= numpoints)) ||
                 (periodic && (index < -1 || index > numpoints)) )
            {
                throw std::out_of_range("Grid index out of range.");
            }
        }

        size_t idx = 0;
        size_t fac = 1;
        for (size_t i=0; i < GridBase<T>::GetDimension(); ++i) {
            if (GridBase<T>::GetPeriodic(i)) {
                idx += indices.at(i) * fac;
                fac *= GridBase<T>::GetNumPoints(i);
            } else {
                idx += (indices.at(i) + 1) * fac;
                fac *= GridBase<T>::GetNumPoints(i) + 2;
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
      : GridBase<T>(numPoints, lower, upper, isPeriodic)
    {
        size_t data_size = 1;
        for (size_t d = 0; d < GridBase<T>::GetDimension(); ++d) {
            size_t storage_size = GridBase<T>::GetNumPoints(d);
            if (!GridBase<T>::GetPeriodic(d)) { storage_size += 2; }
            data_size *= storage_size;
        }

        GridBase<T>::data_.resize(data_size);
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
        return GridBase<T>::data_.begin();
    }

    //! Return iterator after last element of internal storage
    typename std::vector<T>::iterator end()
    {
        return GridBase<T>::data_.end();
    }

    //! Return const iterator at first element of internal storage
    typename std::vector<T>::const_iterator begin() const
    {
        return GridBase<T>::data_.begin();
    }

    //! Return const iterator after last element of internal storage
    typename std::vector<T>::const_iterator end() const
    {
        return GridBase<T>::data_.end();
    }
};

} // End namespace SSAGES
