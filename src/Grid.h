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
    std::vector<size_t> numPoints_;

    //! Edges of the Grid in each dimension.
    std::pair< std::vector<double>,std::vector<double> > edges_;

    //! Periodicity of the Grid.
    std::vector<bool> isPeriodic_;

    //! Internal storage of the data
    std::vector<T> data_;

public:
    //! Constructor
    /*!
     * \param numPoints Number of grid points in each dimension.
     * \param lower Lower edges of the grid.
     * \param upper Upper edges of the grid.
     * \param isPeriodic Bools specifying the periodicity in the respective
     *                   dimension. Default: Non-periodic in all dimensions.
     *
     * The constructor is intentionally private to make sure, that it is only
     * called via BuildGrid().
     *
     * The dimension of the grid is determined by the size of the parameter
     * vectors.
     */
    Grid(std::vector<size_t> numPoints,
         std::vector<double> lower,
         std::vector<double> upper,
         std::vector<bool> isPeriodic = std::vector<bool>())
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
            data_size *= GetNumPoints(d);
        }

        data_.resize(data_size);
    }

    //! Get the dimension.
    size_t GetDimension() const
    {
        return dimension_;
    }

    //! Get the number of points for all dimensions.
    std::vector<size_t> GetNumPoints() const
    {
        return numPoints_;
    }

    //! Get the number of points for a specific dimension.
    /*!
     * \param dim Index of the dimension.
     *
     * \note The first dimension uses the index 0.
     */
    size_t GetNumPoints(size_t dim) const
    {
        if (dim >= GetDimension()) {
            std::cerr << "Warning! Grid size requested for a dimension larger "
                         "than the grid dimensionality!\n";
            return 0;
        }

        return numPoints_.at(dim);
    }

    //! Return the lower edges of the Grid.
    std::vector<double> GetLower() const
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
    std::vector<double> GetUpper() const
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
    std::vector<bool> GetPeriodic() const
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
    std::vector<int> GetIndices(std::vector<double> x) const
    {
        return std::vector<int>(GetDimension(), 0);
    }

    //! Access Grid element read-only
    /*!
     * \param indices Vector of integers specifying the grid point.
     */
    const T& at(std::vector<int> indices) const
    {
        return data_.at(0);
    }

    //! Access Grid element read/write
    /*!
     * \param indices Vector of integers specifying the grid point.
     */
    T& at(std::vector<int> indices)
    {
        return data_.at(0);
    }

    //! Access Grid element pertaining to a specific point -- read-only
    /*!
     * \param x Vector of doubles specifying a point.
     *
     * This function is provided for convenience. It is identical to
     * Grid::at(Grid::GetIndices(x)).
     */
    const T& at(std::vector<double> x) const
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
    T& at(std::vector<double> x)
    {
        return data_.at(GetIndices(x));
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
        std::vector<size_t> number_points;
        for (auto &np : json["number_points"]) {
            number_points.push_back(np.asUInt());
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
        Grid<T>* grid = new Grid(number_points, lower, upper);

        // Since BuildGrid() is a method of Grid, we can access private
        // variables.
        grid->isPeriodic_ = isPeriodic;

        return grid;
    }

    //! \copydoc Serializable::Serialize()
    /*!
     * \warning Serialization not yet implemented.
     */
    void Serialize(Json::Value& json) const override
    {

    }
};

} // End namespace SSAGES
