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
 * @file    sort.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the interface for parallel sorting.
 *
 * TODO:
 * - [x] fix stable sort
 * - [x] fix sort on non GCC compilers
 * - [ ] radix sort
 * - [ ] fix sorting of samples
 * - [ ] implement and try out different parallel sorting algorithms
 * - [ ] bitonic sort for single elemenet per processor
 */

#ifndef MXX_SORT_HPP
#define MXX_SORT_HPP

#include "comm_fwd.hpp"
#include "samplesort.hpp"

namespace mxx {

template<typename _Iterator, typename _Compare>
void sort(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm = mxx::comm()) {
    // use sample sort
    impl::samplesort<_Iterator, _Compare, false>(begin, end, comp, comm);
}

template <typename _Iterator>
void sort(_Iterator begin, _Iterator end, const mxx::comm& comm = mxx::comm()) {
    typedef std::less<typename std::iterator_traits<_Iterator>::value_type> Cmp;
    impl::samplesort<_Iterator, Cmp, false>(begin, end, Cmp(), comm);
}

template<typename _Iterator, typename _Compare>
void stable_sort(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm = mxx::comm()) {
    // use stable sample sort
    impl::samplesort<_Iterator, _Compare, true>(begin, end, comp, comm);
}

template <typename _Iterator>
void stable_sort(_Iterator begin, _Iterator end, const mxx::comm& comm = mxx::comm()) {
    typedef std::less<typename std::iterator_traits<_Iterator>::value_type> Cmp;
    impl::samplesort<_Iterator, Cmp, true>(begin, end, Cmp(), comm);
}

template<typename _Iterator, typename _Compare>
bool is_sorted(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm = mxx::comm()) {
    return impl::is_sorted(begin, end, comp, comm);
}

template<typename _Iterator>
bool is_sorted(_Iterator begin, _Iterator end, const mxx::comm& comm = mxx::comm()) {
    typedef std::less<typename std::iterator_traits<_Iterator>::value_type> Cmp;
    return impl::is_sorted(begin, end, Cmp(), comm);
}

// assumes input is sorted, removes duplicates in global range
template <typename Iterator, typename BinaryPredicate>
Iterator unique(Iterator begin, Iterator end, BinaryPredicate eq, const mxx::comm& comm = mxx::comm()) {
    typedef typename std::iterator_traits<Iterator>::value_type T;
    Iterator dest = begin;
    mxx::comm c = comm.split(begin != end);
    comm.with_subset(begin != end, [&](const mxx::comm& c) {
        size_t n = std::distance(begin, end);

        // send last item to next processor
        T last = *(begin + (n-1));
        T prev = mxx::right_shift(last, c);

        // skip elements which are equal to the last one on the previous processor
        if (c.rank() > 0)
            while (begin != end && eq(prev, *begin))
                ++begin;
        if (begin == end)
            return;
        *dest = *begin;

        // remove duplicates
        while (++begin != end)
            if (!eq(*dest, *begin))
                *++dest = *begin;
        ++dest;
    });
    return dest;
}

template <typename Iterator>
Iterator unique(Iterator begin, Iterator end, const mxx::comm& comm = mxx::comm()) {
    return unique(begin, end, std::equal_to<typename std::iterator_traits<Iterator>::value_type>(), comm);
}

#include "comm_def.hpp"

} // namespace mxx

#endif // MXX_SORT_HPP
