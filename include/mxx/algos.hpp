/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    algos.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @author  Nagakishore Jammula <njammula3@mail.gatech.edu>
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements some common sequential algorithms
 *
 */

#ifndef MXX_ALGOS_H
#define MXX_ALGOS_H

#include <cstdlib>
#include <vector>

#include "common.hpp"

namespace mxx {

/**
 * @brief  Calculates the inclusive prefix sum of the given input range.
 *
 * @param begin An iterator to the beginning of the sequence.
 * @param end   An iterator to the end of the sequence.
 */
template <typename Iterator>
void prefix_sum(Iterator begin, const Iterator end)
{
    // set the total sum to zero
    typename std::iterator_traits<Iterator>::value_type sum = 0;

    // calculate the inclusive prefix sum
    while (begin != end)
    {
        sum += *begin;
        *begin = sum;
        ++begin;
    }
}

/**
 * @brief  Calculates the exclusive prefix sum of the given input range.
 *
 * @param begin An iterator to the beginning of the sequence.
 * @param end   An iterator to the end of the sequence.
 */
template <typename Iterator>
void excl_prefix_sum(Iterator begin, const Iterator end)
{
    // set the total sum to zero
    typename std::iterator_traits<Iterator>::value_type sum = 0;
    typename std::iterator_traits<Iterator>::value_type tmp;

    // calculate the exclusive prefix sum
    while (begin != end)
    {
        tmp = sum;
        sum += *begin;
        *begin = tmp;
        ++begin;
    }
}

#if 0

/**
 * @brief Returns whether the given input range is in sorted order.
 *
 * @param begin An iterator to the begin of the sequence.
 * @param end   An iterator to the end of the sequence.
 *
 * @return  Whether the input sequence is sorted.
 */
template <typename Iterator>
bool is_sorted(Iterator begin, Iterator end)
{
    if (begin == end)
        return true;
    typename std::iterator_traits<Iterator>::value_type last = *(begin++);
    while (begin != end)
    {
        if (last > *begin)
            return false;
        last = *(begin++);
    }
    return true;
}

/**
 * @brief Performs a balanced partitioning.
 *
 * This function internally does 3-way partitioning (i.e. partitions the
 * input sequence into elements of the three classes given by:  <, ==, >
 *
 * After partioning, the elements that are in the `==` class are put into the
 * class of `<` or `>`.
 *
 * This procedure ensures that if the pivot is part of the sequence, then
 * the returned left (`<`) and right (`>`) sequences are both at least one
 * element in length, i.e. neither is of length zero.
 *
 * @param begin An iterator to the beginning of the sequence.
 * @param end   An iterator to the end of the sequence.
 * @param pivot A value to be used as pivot.
 *
 * @return An iterator pointing to the first element of the right (`>`)
 *         sequence.
 */
template<typename Iterator, typename T>
Iterator balanced_partitioning(Iterator begin, Iterator end, T pivot)
{
    if (begin == end)
        return begin;

    // do 3-way partitioning and then balance the result
    Iterator eq_it = begin;
    Iterator le_it = begin;
    Iterator ge_it = end - 1;

    while (true)
    {
        while (le_it < ge_it && *le_it < pivot) ++le_it;
        while (le_it < ge_it && *ge_it > pivot) --ge_it;

        if (le_it == ge_it)
        {
            if (*le_it == pivot)
            {
                std::swap(*(le_it++), *(eq_it++));
            }
            else if (*le_it < pivot)
            {
                le_it++;
            }
            break;
        }

        if (le_it > ge_it)
            break;

        if (*le_it == pivot)
        {
            std::swap(*(le_it++), *(eq_it++));
        } else if (*ge_it == pivot)
        {
            std::swap(*(ge_it--), *(le_it));
            std::swap(*(le_it++), *(eq_it++));
        }
        else
        {
            std::swap(*(le_it++), *(ge_it--));
        }
    }


    // if there are pivots, put them to the smaller side
    if (eq_it != begin)
    {
        typedef typename std::iterator_traits<Iterator>::difference_type diff_t;
        diff_t left_size = le_it - eq_it;
        diff_t right_size = end - le_it;

        // if left is bigger, put to right, otherwise just leave where they are
        if (left_size > right_size)
        {
            // put elements to right
            while (eq_it > begin)
            {
                std::swap(*(--le_it), *(--eq_it));
            }
        }
    }
    // no pivots, therefore simply return the partitioned sequence
    return le_it;
}

#endif

/*********************************************************************
 *                       Bucketing algorithms                        *
 *********************************************************************/

// TODO: iterator version?

/**
 * @brief   Inplace bucketing of values into `num_buckets` buckets.
 *
 * This particular implementation uses an internal temporary buffer
 * of the same size as the input. Thus requiring that amount of additional
 * memory space. For a version that doesn't need O(n) additional memory,
 * use the (somewhat slower) `bucketing_inplace()` function below.
 *
 * @tparam T            Input type
 * @tparam Func         Type of the key function.
 * @param input[in|out] Contains the values to be bucketed.
 *                      This vector is both input and output.
 * @param key_func      A function taking a type T and returning the bucket index
 *                      in the range [0, num_buckets).
 * @param num_buckets   The total number of buckets.
 *
 * @return              The number of elements in each bucket.
 *                      The size of this vector is `num_buckets`.
 */
template <typename T, typename Func>
std::vector<size_t> bucketing(std::vector<T>& input, Func key_func, size_t num_buckets) {
    // initialize number of elements per bucket
    std::vector<size_t> bucket_counts(num_buckets, 0);

    // if input is empty, simply return
    if (input.size() == 0)
        return bucket_counts;

    // [1st pass]: counting the number of elements per bucket
    for (auto it = input.begin(); it != input.end(); ++it) {
        MXX_ASSERT(0 <= key_func(*it) && (size_t)key_func(*it) < num_buckets);
        bucket_counts[key_func(*it)]++;
    }

    // get offsets of where buckets start (= exclusive prefix sum)
    std::vector<std::size_t> offset(bucket_counts.begin(), bucket_counts.end());
    excl_prefix_sum(offset.begin(), offset.end());

    // [2nd pass]: saving elements into correct position
    std::vector<T> tmp_result(input.size());
    for (auto it = input.begin(); it != input.end(); ++it) {
        tmp_result[offset[key_func(*it)]++] = *it;
    }

    // replacing input with tmp result buffer and return the number of elements
    // in each bucket
    input.swap(tmp_result);
    return bucket_counts;
}

// inplace version (doesn't require O(n) additional memory like the other two approaches)
/**
 * @brief   Inplace bucketing of values into `num_buckets` buckets.
 *
 * This particular implementation is truly inplace, and doesn't require O(n)
 * additional memory like the `bucketing()` function above.
 * However, this implementation is slightly slower (~1.5x).
 *
 * @tparam T            Input type
 * @tparam Func         Type of the key function.
 * @param input[in|out] Contains the values to be bucketed.
 *                      This vector is both input and output.
 * @param key_func      A function taking a type T and returning the bucket index
 *                      in the range [0, num_buckets).
 * @param num_buckets   The total number of buckets.
 *
 * @return              The number of elements in each bucket.
 *                      The size of this vector is `num_buckets`.
 */
template <typename T, typename Func>
std::vector<size_t> bucketing_inplace(std::vector<T>& input, Func key_func, size_t num_buckets) {
    // initialize number of elements per bucket
    std::vector<size_t> bucket_counts(num_buckets, 0);

    // if input is empty, simply return
    if (input.size() == 0)
        return bucket_counts;

    // [1st pass]: counting the number of elements per bucket
    for (auto it = input.begin(); it != input.end(); ++it) {
        MXX_ASSERT(0 <= key_func(*it) && (size_t)key_func(*it) < num_buckets);
        bucket_counts[key_func(*it)]++;
    }
    // get exclusive prefix sum
    // get offsets of where buckets start (= exclusive prefix sum)
    // and end (=inclusive prefix sum)
    std::vector<size_t> offset(bucket_counts.begin(), bucket_counts.end());
    std::vector<size_t> upper_bound(bucket_counts.begin(), bucket_counts.end());
    excl_prefix_sum(offset.begin(), offset.end());
    prefix_sum(upper_bound.begin(), upper_bound.end());

    // in-place bucketing
    size_t cur_b = 0;
    for (size_t i = 0; i < input.size();) {
        // skip full buckets
        while (cur_b < num_buckets-1 && offset[cur_b] >= upper_bound[cur_b]) {
            // skip over full buckets
            i = offset[++cur_b];
        }
        // break if all buckets are done
        if (cur_b >= num_buckets-1)
            break;
        size_t target_b = key_func(input[i]);
        MXX_ASSERT(0 <= target_b && target_b < num_buckets);
        if (target_b == cur_b) {
            // item correctly placed
            ++i;
        } else {
            // swap to correct bucket
            MXX_ASSERT(target_b > cur_b);
            std::swap(input[i], input[offset[target_b]]);
        }
        offset[target_b]++;
    }
    return bucket_counts;
}

// Tony's old implementation for bucketing:
// (needs more memory and is slower than the two versions above)

// Complexity is O(b)*O(n), where as the two version above are O(b + n)
// scales badly, same as bucketing_copy, but with a factor of 2 for large
// data.
// when fixed n, slight increase with b..  else increases with n.
template<typename T, typename Func>
std::vector<size_t> bucketing_tony(std::vector<T>& elements, Func key_func, size_t num_buckets) {

    // number of elements per bucket
    std::vector<size_t> send_counts(num_buckets, 0);

    // if no elements, return 0 count for each bucket
    if (elements.size() == 0)
        return send_counts;

    // for each element, track which bucket it belongs into
    std::vector<long> membership(elements.size());
    for (size_t i = 0; i < elements.size(); ++i)
    {
        membership[i] = key_func(elements[i]);
        ++(send_counts[membership[i]]);
    }
    // at this point, have target assignment for each data element, and also count for each process bucket.

    // compute the offsets within the buffer
    std::vector<size_t> offset = send_counts;
    excl_prefix_sum(offset.begin(), offset.end());
    std::vector<size_t> maxes = offset;

    for (size_t i = 0; i < num_buckets; ++i) {
        maxes[i] += send_counts[i];
    }


    //== swap elements around.
    T val;
    size_t tar_pos, start_pos;

    long target;

    // while loop will stop under 2 conditions:
    //      1. returned to starting position (looped), or
    //      2, tar_pos is the current pos.
    // either way, we need a new starting point.  instead of searching through buffer O(N), search
    // for incomplete buckets via offset O(p).

    for (size_t i = 0; i < num_buckets;) {
        // determine the starting position.
        if (offset[i] == maxes[i]) {
            ++i;  // skip all completed buckets
            continue;  // have the loop check value.
        }
        // get the start pos.
        start_pos = offset[i];

        // set up the variable with the current entry.
        target = membership[start_pos];
        if (target > -1) {
            val = ::std::move(elements[start_pos]);  // value to move
            membership[start_pos] = -2;                // special value to indicate where we started from.

            while (target > -1) {  // if -1 or -2, then either visited or beginning of chain.
                tar_pos = offset[target]++;  // compute new position.  earlier offset values for the same pid are should have final values already.
                target = membership[tar_pos];

                // save the info at tar_pos;
                ::std::swap(val, elements[tar_pos]);  // put what's in src into buffer at tar_pos, and save what's at buffer[tar_pos]
                membership[tar_pos] = -1;               // mark as visited.

            }  // else already visited, so done.
        }
    }

    return send_counts;
}

} // namespace mxx

#endif // MXX_ALGOS_H
