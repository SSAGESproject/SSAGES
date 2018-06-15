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
#include <fstream>
#include <istream>

#include "Drivers/DriverException.h"
#include "GridBase.h"
#include "schema.h"
#include "Validator/ObjectRequirement.h"

namespace SSAGES
{

//! Basic Grid.
/*!
 * \tparam T type of data to be stored in the grid.
 *
 * A Grid is a method to store data in SSAGES. It is used to discretize a
 * continuous number, typically a collective variable or a function, into
 * \c number_points grid points. At each grid point, an arbitrary type of
 * data can be stored, specified via the template parameter \c T.
 *
 * The grid can be of arbitrary dimension. For each dimension, the lower
 * bound, the upper bound and the number of grid points need to be specified.
 * Furthermore, the grid can be defined as periodic or non-periodic in the
 * respective dimension. By default, the grid is non-periodic. The grid points
 * are indexed from 0 to number_points-1 following the standard C/C++
 * convention.
 *
 * The grid spacing \c Delta is given by (upper - lower)/number_points. Thus,
 * grid point \c n corresponds to the interval [lower + n*Delta, lower + (n+1)*Delta).
 * Note that n follows the C/C++ convention, i.e. n = 0 for the first interval.
 * The bin indices pertaining to a given point can be obtained via GetIndices().
 *
 * \ingroup Core
 */
template<typename T>
class Grid : public GridBase<T>
{
private:
    //! Map d-dimensional indices to 1-d data vector
    /*!
     * \param indices Vector specifying the grid point.
     * \return Index of 1d data vector.
     *
     * Map a set of indices to the index of the 1d data vector. Keep in mind,
     * that the data includes underflow (index -1) and overflow (index
     * numPoints) bins in periodic dimension.
     */
    size_t mapTo1d(const std::vector<int> &indices) const override
    {
        // Check if an index is out of bounds
        for (size_t i = 0; i < GridBase<T>::GetDimension(); ++i) {
            int index = indices.at(i);
            int numpoints = GridBase<T>::numPoints_[i];
            if ( index < 0 || index >= numpoints )
            {
                throw std::out_of_range("Grid index out of range.");
            }
        }

        size_t idx = 0;
        size_t fac = 1;
        for (size_t i = 0; i < GridBase<T>::GetDimension(); ++i) {
            idx += indices.at(i) * fac;
            fac *= GridBase<T>::numPoints_[i];
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
            size_t storage_size = GridBase<T>::numPoints_[d];
            data_size *= storage_size;
        }

        GridBase<T>::data_.resize(data_size);
    }

    //! Set up the grid
    /*!
     * \param json JSON value containing all input information.
     * \return Pointer to the newly built grid.
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
     * \return Pointer to the newly built grid.
     *
     * This function builds a grid from a JSON node. It will return a nullptr
     * if an unknown error occured, but generally, it will throw a
     * BuildException on failure.
     */
    static Grid<T>* BuildGrid(const Json::Value& json, const std::string& path)
    {
		Json::ObjectRequirement validator;
		Json::Value schema;
		Json::CharReaderBuilder rbuilder;
		Json::CharReader* reader = rbuilder.newCharReader();

		reader->parse(JsonSchema::grid.c_str(),
		              JsonSchema::grid.c_str() + JsonSchema::grid.size(),
		              &schema, NULL);
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

    //! Custom Iterator
    /*!
     * This iterator is designed of travesing through a grid. The starting point
     * is at grid index 0 for each dimension. The last valid grid point has the
     * num_points in each dimension, where num_points is the number of grid
     * points in the respective dimension.
     *
     * The iterator can be used as a standard iterator with operator* accessing
     * the grid point at which the iterator currently is.
     *
     * Additionally, the functions GridIterator::indices() and
     * GridIterator::coordinates() are provided. These functions return the
     * indices of the current grid point and the center of the grid point
     * interval in real space, respectively.
     *
     * The iterator can be moved to an arbitrary position. As indices() returns
     * a reference (and not a const reference), it can be used to move the
     * iterator. For example:
     *
     * \code{.cpp}
     * HistIterator it = hist->begin();
     * it.indices() = {1,1,1};
     * \endcode
     *
     * moves the iterator to the grid point [1, 1, 1].
     *
     * The iterator can be traversed in a standard fashion with the increment
     * and decrement operators operator++ and operator--. When the increment
     * operator is invoked, the bin index for the lowest dimension is increased
     * by 1. If it moves beyond the valid range in this dimension, the index is
     * reset to 0 and the index of the next higher dimension is increased by 1.
     * The decrement operator traveses the grid in the same fashion but opposite
     * direction.
     *
     * Additionaly, the iterator can be shifted by adding or subtracting a vector
     * of ints. The vector needs to have the same dimension as the histogram.
     */
    template<typename R>
    class GridIterator {
    public:
        //! Type name of the iterator.
        typedef GridIterator self_type;

        //! Difference type is an int.
        typedef int difference_type;

        //! Either T or const T for iterator and const_iterator, respectively.
        typedef R value_type;

        //! Either T* or T const* for iterator and const_iterator, respectively.
        typedef R* pointer;

        //! Either T& or T const& for iterator and const_iterator, respectively.
        typedef R& reference;

        //! HistIterator is a bidirectional iterator.
        typedef std::bidirectional_iterator_tag iterator_category;

        //! Use default constructor.
        GridIterator() = default;

        //! Constructor
        /*!
         * \param indices Bin indices specifying the current position of the
         *                iterator.
         * \param grid Pointer to the grid to iterate over.
         */
        GridIterator(const std::vector<int> &indices, Grid<T> *grid)
            : indices_(indices), grid_(grid)
        {
        }

        //! Const constructor
        /*!
         * \param indices Bin indices specifying the current position of the
         *                iterator.
         * \param grid Pointer to the grid to iterate over.
         */
        GridIterator(const std::vector<int> &indices, const Grid<T> *grid)
            : indices_(indices), grid_(grid)
        {
        }

        //! Copy constructor
        /*!
         * \param other GridIterator to be copied.
         */
        GridIterator(const GridIterator &other)
            : indices_(other.indices_), grid_(other.grid_)
        {
        }

        //! Dereference operator.
        /*!
         * \return Reference to the value at the current grid position.
         */
        reference operator*() { return grid_->at(indices_); }

        //! Pre-increment operator.
        /*!
         * \return Reference to iterator.
         *
         * Increments the bin index of lowest dimension. If an index moves
         * beyond the maximum value (num_points-1), it is reset to 0 and the
         * index of the next higher dimension is increased by 1.
         */
        self_type &operator++()
        {
            indices_.at(0) += 1;
            for (size_t i = 0; i < grid_->GetDimension() - 1; ++i) {
                if (indices_.at(i) >= grid_->numPoints_[i]) {
                    indices_.at(i) = 0;
                    indices_.at(i+1) += 1;
                }
            }

            return *this;
        }

        //! Post-increment operator.
        /*!
         * \return Copy of iterator before incrementing.
         */
        self_type operator++(int)
        {
            GridIterator it(*this);
            ++(*this);
            return it;
        }

        //! Addition assignment operator
        /*
         * \param shift Vector of shifts in each dimension.
         * \return Reference to itself.
         *
         * This operator shifts the current position of the iterator by the
         * given amount in each dimension.
         *
         * Example:
         *
         * \code{.cpp}
         * if += {1,1,1};
         * \endcode
         *
         * In this example the current position of the iterator is shifted
         * diagonally.
         */
        self_type &operator+=(std::vector<int> shift)
        {
            if (shift.size() != grid_->GetDimension()) {
                throw std::invalid_argument("Vector to shift iterator does not "
                                            "match grid dimension.");
            }

            for (size_t i = 0; i < grid_->GetDimension(); ++i) {
                indices_.at(i) += shift.at(i);
            }

            return *this;
        }

        //! Addition operator.
        /*!
         * \param shift Amount of shift in each dimension.
         * \return Copy of iterator after shift.
         *
         * Shift the iterator by a given vector.
         */
        const self_type operator+(std::vector<int> shift)
        {
            return GridIterator(*this) += shift;
        }

        //! Pre-decrement operator.
        /*!
         * \return Reference to iterator after decrementing.
         *
         * Traveses the histogram in the opposite direction to the increment
         * operator.
         */
        self_type &operator--()
        {
            indices_.at(0) -= 1;
            for (size_t i = 0; i < grid_->GetDimension() - 1; ++i) {
                if (indices_.at(i) < 0) {
                    indices_.at(i) = grid_->numPoints_[i]-1;
                    indices_.at(i+1) -= 1;
                }
            }

            return *this;
        }

        //! Post-decrement operator.
        /*!
         * \return Copy of iterator before decrementing.
         */
        self_type operator--(int)
        {
            GridIterator it(*this);
            --(*this);
            return it;
        }

        //! Subtraction assignment operator.
        /*!
         * \param shift Vector to be subtracted from the current grid indices.
         * \return Reference to iterator.
         */
        self_type &operator-=(std::vector<int> shift)
        {
            if (shift.size() != grid_->GetDimension()) {
                throw std::invalid_argument("Vector to shift iterator does not "
                                            "match histogram dimension.");
            }

            for (size_t i = 0; i < grid_->GetDimension(); ++i) {
                indices_.at(i) -= shift.at(i);
            }

            return *this;
        }

        //! Subtraction iterator
        /*!
         * \param shift Vector to be subtracted from the current grid indices.
         * \return Copy of iterator after shift.
         */
        const self_type operator-(std::vector<int> shift)
        {
            return GridIterator(*this) -= shift;
        }

        //! Equality operator
        /*!
         * \param rhs Iterator to which this iterator is compared.
         * \return \c True if both iterators access the same grid point on the
         *         same grid. Else return \c False.
         */
        bool operator==(const self_type &rhs) const
        {
            return indices_ == rhs.indices_ && grid_ == rhs.grid_;
        }

        //! Non-equality operator.
        /*!
         * \param rhs Iterator to which this iterator is compared.
         * \return \c False if both iterators access the same grid point on the
         *         same grid. Else return \c True.
         */
        bool operator!=(const self_type &rhs) const
        {
            return !( (*this) == rhs );
        }

        //! Access indices.
        /*!
         * \return Indices of current grid point.
         *
         * \note This function returns a reference and can be used to move the
         *       current grid point.
         */
        std::vector<int> &indices()
        {
            return indices_;
        }

        //! Access a specific index.
        /*!
         * \param d Dimension of the index.
         * \return Index of the current grid point in the specified dimension.
         *
         */
        int &index(size_t d)
        {
            return indices()[d];
        }

        //! Access coordinates.
        /*!
         * \return Center point of the current grid point.
         */
        std::vector<double> coordinates() const
        {
            return grid_->GetCoordinates(indices_);
        }

        //! Access specific coordinate dimension.
        /*!
         * \param d Dimension of the coordinate.
         * \return Center of the current grid point in the specified dimension.
         */
        double coordinate(size_t d) const
        {
            return coordinates()[d];
        }

    private:
        //! Indices of current grid point.
        std::vector<int> indices_;

        //! Pointer to grid to iterate over.
        Grid<T> *grid_;
    };

    //! Custom iterator over a grid.
    typedef GridIterator<T> iterator;

    //! Custom constant iterator over a grid.
    typedef GridIterator<const T> const_iterator;

    //! Return iterator at first grid point.
    /*!
     * \return Iterator at first grid point.
     *
     * The first grid point is defined as the grid point with index 0 in all
     * dimensions.
     */
    iterator begin()
    {
        std::vector<int> indices(GridBase<T>::GetDimension(), 0);

        return iterator(indices, this);
    }

    //! Return iterator after last valid grid point.
    /*!
     * \return Iterator after last valid grid point.
     *
     * The last valid grid point has index num_points - 1 in all dimensions.
     */
    iterator end()
    {
        std::vector<int> indices(GridBase<T>::GetDimension());
        for (size_t i = 0; i < indices.size(); ++i) {
            indices.at(i) = GridBase<T>::numPoints_[i] - 1;
        }

        iterator it(indices, this);
        return ++it;
    }

    //! Return const iterator at first grid point.
    /*!
     * \return Const iterator at first grid point.
     *
     * The first grid point is defined as the grid point with index 0 in all
     * dimensions.
     */
    const_iterator begin() const
    {
        std::vector<int> indices(GridBase<T>::GetDimension(), 0);

        return iterator(indices, this);
    }

    //! Return const iterator after last valid grid point.
    /*!
     * \return Const iterator after last valid grid point.
     *
     * The last valid grid point has index num_points - 1 in all dimensions.
     */
    const_iterator end() const
    {
        std::vector<int> indices(GridBase<T>::GetDimension());
        for (size_t i = 0; i < indices.size(); ++i) {
            indices.at(i) = GridBase<T>::numPoints_[i] - 1;
        }

        iterator it(indices, this);
        return ++it;
    }

    //! Write grid out to file.
    /*!
     * \param filename Name of the output file. 
     */
    void WriteToFile(const std::string& filename)
    {
        std::ofstream output(filename.c_str(), std::ofstream::out);

        // Print out header.
        output << "#! type grid\n";
        output << "#! dim  " << GridBase<T>::dimension_ << "\n"; 
        output << "#! count "; 
        for(auto& c : this->GetNumPoints()) 
            output << c << " "; 
        output << "\n";
        output << "#! lower "; 
        for(auto& l : this->GetLower())
            output << l << " "; 
        output << "\n"; 
        output << "#! upper "; 
        for(auto& u : this->GetUpper())
            output << u << " "; 
        output << "\n"; 
        output << "#! periodic "; 
        for(auto p : GridBase<T>::isPeriodic_)
            output << p << " "; 
        output << "\n";

        for(auto it = this->begin(); it != this->end(); ++it)
        {
            auto coords = it.coordinates(); 
            for(auto& c : coords)
                output << std::setprecision(8) << std::fixed << c << " "; 
            output.unsetf(std::ios_base::fixed);
            output << std::setprecision(16) << *it << "\n";
        }

        output.close();
    }

    //! Builds a grid from file.
    /*!
     * \param filename Filename containing grid data.
     * \return Pointer to the newly built grid.
     *
     * This function builds a grid from a text data file. 
     */
    void LoadFromFile(const std::string& filename)
    {
        std::ifstream file(filename);
        std::string line, buff; 

        // Skip type for now. 
        std::getline(file, line); 

        // Get dimension.
        {
            size_t dim = 0;
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> buff >> buff >> dim; 
            if(dim != GridBase<T>::dimension_)
                throw std::invalid_argument(filename + 
                    ": Expected dimension " + std::to_string(GridBase<T>::dimension_) + 
                    " but got " + std::to_string(dim) + " instead.");
        }

        // Get size. 
        {
            std::getline(file, line); 
            std::istringstream iss(line);
            int count = 0;
            auto counts = GridBase<T>::numPoints_;
            iss >> buff >> buff;
            for(int i = 0; iss >> count; ++i)
            {
				int period = 0;
				if(GridBase<T>::isPeriodic_[i])
				{
					period++;
				}				
                if((count+period) != counts[i])
                    throw std::invalid_argument(filename + 
                    ": Expected count " + std::to_string(counts[i]) + 
                    " on dimension " + std::to_string(i) + " but got " + 
                    std::to_string(count) + " instead.");
            }
        }

        // Get lower bounds. 
        {
            std::getline(file, line); 
            std::istringstream iss(line);
            double lower = 0; 
            auto lowers = GridBase<T>::edges_.first; 
            iss >> buff >> buff;
            for(int i = 0; iss >> lower; ++i)
            {
				double spacing = 0;
				if(GridBase<T>::isPeriodic_[i])
				{
					spacing = (GridBase<T>::edges_.second[i] - GridBase<T>::edges_.first[i]) / GridBase<T>::numPoints_[i];
				}				
                if(std::abs(lower-spacing/2 - lowers[i]) > 1e-8)
                    throw std::invalid_argument(filename + 
                    ": Expected lower " + std::to_string(lowers[i]+spacing/2) + 
                    " on dimension " + std::to_string(i) + " but got " + 
                    std::to_string(lower) + " instead.");
            }
        }

        // Get upper bounds. 
        {
            std::getline(file, line); 
            std::istringstream iss(line);
            double upper = 0; 
            auto uppers = GridBase<T>::edges_.second; 
            iss >> buff >> buff;
            for(int i = 0; iss >> upper; ++i)
            {
				double spacing = 0;
				if(GridBase<T>::isPeriodic_[i])
				{
					spacing = (GridBase<T>::edges_.second[i] - GridBase<T>::edges_.first[i]) / GridBase<T>::numPoints_[i];
				}
                if(std::abs((upper+spacing/2) - uppers[i]) > 1e-8)
                    throw std::invalid_argument(filename + 
                    ": Expected upper " + std::to_string(uppers[i]-spacing/2) + 
                    " on dimension " + std::to_string(i) + " but got " + 
                    std::to_string(upper) + " instead.");
            }
        }

        // Get periodic. 
        {
            std::getline(file, line); 
            std::istringstream iss(line);
            bool periodic; 
            auto periodics = GridBase<T>::isPeriodic_; 
            iss >> buff >> buff;
            for(int i = 0; iss >> periodic; ++i)
            {
                if(periodic != periodics[i])
                    throw std::invalid_argument(filename + 
                    ": Expected periodic " + std::to_string(periodics[i]) + 
                    " on dimension " + std::to_string(i) + " but got " + 
                    std::to_string(periodic) + " instead.");
            }
        }

        // Finally load data.
        for(auto& d : GridBase<T>::data_)
        {
            std::getline(file, line);
            std::istringstream iss(line);
            
            // Skip coordinates.
            for(size_t i = 0; i < GridBase<T>::dimension_; ++i) 
                iss >> buff;
            
            iss >> d;
        }
    }
};

} // End namespace SSAGES
