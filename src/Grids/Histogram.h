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
#include "schema.h"
#include "Validator/ObjectRequirement.h"

namespace SSAGES
{

//! Basic Histogram.
/*!
 * \tparam T type of data to be stored in the histogram.
 *
 * A Histogram is a method to store data in SSAGES. It is used to discretize
 * a continuous number, typically a collective variable, into \c number_points
 * bins. For each bin, an arbitrary type of data can be stored, specified via
 * the template parameter \c T.
 *
 * The histogram can be of arbitrary dimension. For each dimension, the lower
 * bound, the upper bound and the number of grid points need to be specified.
 * Furthermore, the histogram can be defined as periodic or non-periodic in the
 * respective dimension. By default, the histogram is non-periodic. The bins
 * are indexed from 0 to number_points-1 following the standard C/C++
 * convention. In contrast to the Grid, a Histogram additionally includes
 * an under- and an overflow bin in each non-periodic dimension. These can be
 * accessed via the indices -1 an \c number_points (see below).
 *
 * The bin width \c Delta is given by (upper - lower)/number_points. Thus, bin
 * \c n corresponds to the interval [lower + n*Delta, lower + (n+1)*Delta).
 * Note that n follows the C/C++ convention, i.e. n = 0 for the first interval.
 * The bin indices pertaining to a given point can be obtained via GetIndices().
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
class Histogram : public GridBase<T>
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
        for (size_t i=0; i < GridBase<T>::GetDimension(); ++i) {
            int index = indices.at(i);
            bool periodic = GridBase<T>::GetPeriodic(i);
            int numpoints = GridBase<T>::GetNumPoints(i);
            if ( (periodic && (index < 0 || index >= numpoints)) ||
                 (periodic && (index < -1 || index > numpoints)) )
            {
                throw std::out_of_range("Bin index out of range.");
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
    Histogram(std::vector<int> numPoints,
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

    //! Set up the histogram
    /*!
     * \param json JSON value containing all input information.
     * \return Pointer to the newly built histogram.
     *
     * This function builds a histogram from a JSON node. It will return a nullptr
     * if an unknown error occured, but generally, it will throw a
     * BuildException of failure.
     */
    static Histogram<T>* BuildHistogram(const Json::Value& json)
    {
        return BuildHistogram(json, "#/Histogram");
    }

    //! Set up the histogram
    /*!
     * \param json JSON Value containing all input information.
     * \param path Path for JSON path specification.
     * \return Pointer to the newly built histogram.
     *
     * This function builds a histogram from a JSON node. It will return a nullptr
     * if an unknown error occured, but generally, it will throw a
     * BuildException on failure.
     */
    static Histogram<T>* BuildHistogram(const Json::Value& json, const std::string& path)
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
        Histogram<T>* hist = new Histogram(number_points, lower, upper, isPeriodic);

        return hist;
    }

    //! Custom Iterator
    /*!
     * This iterator is designed of travesing through a histogram. The starting
     * point is at grid index 0 (-1) for each periodic (non-periodic) dimension.
     *
     * The iterator can be used as a standard iterator with operator* accessing
     * the grid point at which the iterator currently is.
     *
     * Additionally, the functions HistIterator::indices() and
     * HistIterator::coordinates() are provided. These functions return the
     * indices of the current bin and the bin center in real space, respectively.
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
     * moves the iterator to the bin with indices [1, 1, 1].
     *
     * The iterator can be traversed in a standard fashion with the increment
     * and decrement operators operator++ and operator--. When the increment
     * operator is invoked, the bin index for the lowest dimension is increased
     * by 1. If it is beyond the histogram size in this dimension, the index is
     * reset to its smallest value (0 for periodic, -1 for non-periodic
     * dimensions) and the index of the next higher dimension is increased by 1.
     * The decrement operator traveses the grid in the same fashion but opposite
     * direction.
     *
     * Additionaly, the iterator can be shifted by adding or subtracting a vector
     * of ints. The vector needs to have the same dimension as the histogram.
     */
    template<typename R>
    class HistIterator {
    public:
        //! Type name of the iterator.
        typedef HistIterator self_type;

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
        HistIterator() = default;

        //! Constructor
        /*!
         * \param indices Bin indices specifying the current position of the
         *                iterator.
         * \param hist Pointer to the histogram to iterate over.
         */
        HistIterator(const std::vector<int> &indices, Histogram<T> *hist)
            : indices_(indices), hist_(hist)
        {
        }

        //! Copy constructor
        /*!
         * \param other HistIterator to be copied.
         */
        HistIterator(const HistIterator &other)
            : indices_(other.indices_), hist_(other.hist_)
        {
        }

        //! Dereference operator.
        /*!
         * \return Reference to the value at the current grid position.
         */
        reference operator*() { return hist_->at(indices_); }

        //! Pre-increment operator.
        /*!
         * \return Reference to iterator.
         *
         * Increments the bin index of the lowest dimension. If an index moves
         * beyond the maximum value (num_points-1 for periodic and num_points
         * for non-periodic dimensions), it is reset to its smallest value (0
         * for periodic, -1 for non-periodic dimensions) and the index of the
         * next higher dimension is increased by 1.
         */
        self_type &operator++()
        {
            indices_.at(0) += 1;
            for (size_t i = 0; i < hist_->GetDimension() - 1; ++i) {
                if (hist_->GetPeriodic(i) &&
                    indices_.at(i) >= hist_->GetNumPoints(i)) {

                    indices_.at(i) = 0;
                    indices_.at(i+1) += 1;
                } else if (!hist_->GetPeriodic(i) &&
                           indices_.at(i) > hist_->GetNumPoints(i)) {

                    indices_.at(i) = -1;
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
            HistIterator it(*this);
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
         * if += {1,1,1};
         *
         * In this example the current position of the iterator is shifted
         * diagonally.
         */
        self_type &operator+=(std::vector<int> shift)
        {
            if (shift.size() != hist_->GetDimension()) {
                throw std::invalid_argument("Vector to shift iterator does not "
                                            "match histogram dimension.");
            }

            for (size_t i = 0; i < hist_->GetDimension(); ++i) {
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
            return HistIterator(*this) += shift;
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
            for (size_t i = 0; i < hist_->GetDimension() - 1; ++i) {
                if (hist_->GetPeriodic(i) && indices_.at(i) < 0) {
                    indices_.at(i) = hist_->GetNumPoints(i)-1;
                    indices_.at(i+1) -= 1;
                } else if (!hist_->GetPeriodic(i) && indices_.at(i) < -1) {
                    indices_.at(i) = hist_->GetNumPoints(i);
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
            HistIterator it(*this);
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
            if (shift.size() != hist_->GetDimension()) {
                throw std::invalid_argument("Vector to shift iterator does not "
                                            "match histogram dimension.");
            }

            for (size_t i = 0; i < hist_->GetDimension(); ++i) {
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
            return HistIterator(*this) -= shift;
        }

        //! Equality operator
        /*!
         * \param rhs Iterator to which this iterator is compared.
         * \return \c True if both iterators access the same grid point on the
         *         same grid. Else return \c False.
         */
        bool operator==(const self_type &rhs) const
        {
            return indices_ == rhs.indices_ && hist_ == rhs.hist_;
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
         * \return Indices of current bin.
         *
         * \note This function returns a reference and can be used to move the
         *       current bin.
         */
        std::vector<int> &indices()
        {
            return indices_;
        }

        //! Access a specific index.
        /*!
         * \param d Dimension of the index.
         * \return Index of the current bin in the specified dimension.
         *
         */
        int &index(size_t d)
        {
            return indices()[d];
        }

        //! Check if current iterator position is under- or overflow bin.
        /*!
         * \return \c True if current bin is an underflow or an overflow bin.
         */
        bool isUnderOverflowBin() const
        {
            for (size_t i = 0; i < indices_.size(); ++i) {
                if (indices_.at(i) == -1 || indices_.at(i) == hist_->GetNumPoints(i)) {
                    return true;
                }
            }
            return false;
        }

        //! Access coordinates.
        /*!
         * \return Center point of the current bin.
         */
        std::vector<double> coordinates() const
        {
            return hist_->GetCoordinates(indices_);
        }

        //! Access specific coordinate dimension.
        /*!
         * \param d Dimension of the coordinate.
         * \return Center point of the current bin in the specified dimension.
         */
        double coordinate(size_t d) const
        {
            return coordinates()[d];
        }

    private:
        //! Indices of current bin.
        std::vector<int> indices_;

        //! Pointer to histogram to iterate over.
        Histogram<T> *hist_;
    };

    //! Custom iterator over a histogram.
    typedef HistIterator<T> iterator;

    //! Custom constant iterator over a histogram.
    typedef HistIterator<const T> const_iterator;

    //! Return iterator at first bin of histogram
    /*!
     * \return Iterator at first bin of the histogram.
     *
     * The first bin is the bin that has the lowest allowed index in all
     * dimensions, i.e. 0 in periodic and -1 in non-periodic dimensions.
     */
    iterator begin()
    {
        std::vector<int> indices(GridBase<T>::GetDimension());
        for (size_t i = 0; i < indices.size(); ++i) {
            if(GridBase<T>::GetPeriodic(i)) {
                indices.at(i) = 0;
            } else {
                indices.at(i) = -1;
            }
        }

        return iterator(indices, this);
    }

    //! Return iterator after last valid bin.
    /*!
     * \return Iterator after last valid bin.
     *
     * The last valid bin is the bin that has the highest allowed index in all
     * dimensions, i.e. num_points - 1 in periodic and num_points in
     * non-periodic dimensions.
     */
    iterator end()
    {
        std::vector<int> indices(GridBase<T>::GetDimension());
        for (size_t i = 0; i < indices.size(); ++i) {
            if (GridBase<T>::GetPeriodic(i)) {
                indices.at(i) = GridBase<T>::GetNumPoints(i) - 1;
            } else {
                indices.at(i) = GridBase<T>::GetNumPoints(i);
            }
        }

        iterator it(indices, this);
        return ++it;
    }

    //! Return const iterator at first bin of histogram
    /*!
     * \return Const iterator at first bin of the histogram.
     *
     * The first bin is the bin that has the lowest allowed index in all
     * dimensions, i.e. 0 in periodic and -1 in non-periodic dimensions.
     */
    typename std::vector<T>::const_iterator begin() const
    {
        std::vector<int> indices(GridBase<T>::GetDimension());
        for (size_t i = 0; i < indices.size(); ++i) {
            if(GridBase<T>::GetPeriodic(i)) {
                indices.at(i) = 0;
            } else {
                indices.at(i) = -1;
            }
        }

        return iterator(indices, this);
    }

    //! Return const iterator after last valid bin.
    /*!
     * \return Const iterator after last valid bin.
     *
     * The last valid bin is the bin that has the highest allowed index in all
     * dimensions, i.e. num_points - 1 in periodic and num_points in
     * non-periodic dimensions.
     */
    typename std::vector<T>::const_iterator end() const
    {
        std::vector<int> indices(GridBase<T>::GetDimension());
        for (size_t i = 0; i < indices.size(); ++i) {
            if (GridBase<T>::GetPeriodic(i)) {
                indices.at(i) = GridBase<T>::GetNumPoints(i) - 1;
            } else {
                indices.at(i) = GridBase<T>::GetNumPoints(i);
            }
        }

        iterator it(indices, this);
        return ++it;
    }
};

} // End namespace SSAGES
