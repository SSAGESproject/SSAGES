/*
 * Copyright 2016 Georgia Institute of Technology
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
 * @file    bitonicsort.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements parallel, MPI sample sort
 */

#ifndef MXX_BITONICSORT_HPP
#define MXX_BITONICSORT_HPP

#include "comm.hpp"
#include "datatypes.hpp"
#include "collective.hpp"

#include <iterator>
#include <algorithm>
#include <vector>
#include <tgmath.h>


namespace mxx {
namespace impl {

constexpr bool asc = false;
constexpr bool desc = true;

template <typename _Iterator, typename _Compare>
void bitonic_split(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm, int partner, int dir) {
    typedef typename std::iterator_traits<_Iterator>::value_type T;
    size_t np = std::distance(begin, end);
    //std::cout << "[Rank " << comm.rank() << "]: split with partner= " << partner << ", dir=" << dir << std::endl;
    std::vector<T> recv_buf(np);
    std::vector<T> merge_buf(np);

    // send/recv
    mxx::datatype dt = mxx::get_datatype<T>();
    MPI_Sendrecv(&*begin, np, dt.type(), partner, 0,
                 &recv_buf[0], np, dt.type(), partner, 0, comm, MPI_STATUS_IGNORE);

    // merge in `dir` direction into merge buffer
    if ((dir == asc && partner > comm.rank()) || (dir == desc && partner < comm.rank())) {
        // forward merge and keep the `np` smaller elements
        _Iterator l = begin;
        typename std::vector<T>::iterator r = recv_buf.begin();
        typename std::vector<T>::iterator o = merge_buf.begin();

        for (size_t i = 0; i < np; ++i) {
            if (comp(*l, *r)) {
                *o = *l;
                ++o;
                ++l;
            } else {
                *o = *r;
                ++o;
                ++r;
            }
        }
    } else {
        // backward merge and keep the `np` larger elements
        _Iterator l = begin+np-1;
        typename std::vector<T>::iterator r = recv_buf.begin()+np-1;
        typename std::vector<T>::iterator o = merge_buf.begin()+np-1;

        for (size_t i = 0; i < np; ++i) {
            if (comp(*l, *r)) {
                *o = *r;
                --o;
                --r;
            } else {
                *o = *l;
                --o;
                --l;
            }
        }
    }
    // copy results into local buffer
    std::copy(merge_buf.begin(), merge_buf.end(), begin);
}

template <typename _Iterator, typename _Compare>
void bitonic_merge(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm, int pbeg, int pend, int dir) {
    MXX_ASSERT(pbeg <= comm.rank() && comm.rank() < pend);

    // get size and terminate at recursive base-case
    int size = pend - pbeg;
    if (size <= 1)
        return;

    // get next greater power of 2
    int p2 = pow(2, ceil(log(size)/log(2)));
    // merge with splits as is done in the power of 2 case
    int pmid = pbeg + p2/2;

    if (comm.rank() < pmid && comm.rank() + p2/2 < pend) {
        // this processor has a partner in the second half
        int partner_rank = comm.rank() + p2/2;
        bitonic_split(begin, end, comp, comm, partner_rank, dir);
        bitonic_merge(begin, end, comp, comm, pbeg, pmid, dir);
    } else if (comm.rank() < pmid) {
        // this process doesn't have a partner but has to recursively
        // participate in the next merge
        bitonic_merge(begin, end, comp, comm, pbeg, pmid, dir);
    } else { // if (comm.rank() >= pmid)
        int partner_rank = comm.rank() - p2/2;
        bitonic_split(begin, end, comp, comm, partner_rank, dir);
        bitonic_merge(begin, end, comp, comm, pmid, pend, dir);
    }
}

template <typename _Iterator, typename _Compare>
void bitonic_sort_rec(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm, int pbeg, int size, bool dir) {

    // get next power of two
    int p2 = pow(2, ceil(log(size)/log(2)));

    // determine where the two sub-recursions are
    int pmid = pbeg + p2/2;
    int pend = pbeg + p2;
    if (pend > comm.size())
        pend = comm.size();
    // recursive base-case
    if (pend - pbeg <= 1)
        return;

    // sort the two subsequences, where the first is always a power of 2
    if (comm.rank() < pmid) {
        // sort descending
        bitonic_sort_rec(begin, end, comp, comm, pbeg, p2/2, !dir);
    } else {
        // sort ascending
        bitonic_sort_rec(begin, end, comp, comm, pmid, size - p2/2, dir);
    }
    // merge bitonic decreasing sequence
    bitonic_merge(begin, end, comp, comm, pbeg, pend, dir);
}

} // namespace impl


template <typename _Iterator, typename _Compare>
void bitonic_sort(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm) {
    size_t np = std::distance(begin, end);
    if (!mxx::all_same(np, comm)) {
        throw std::runtime_error("bitonic sort only works for the same number of elements on each process.");
    }

    if (!std::is_sorted(begin, end, comp)) {
        // sort local elements first
        std::sort(begin, end, comp);
    }

    // recursively sort via bitonic sort
    impl::bitonic_sort_rec(begin, end, comp, comm, 0, comm.size(), impl::asc);
}

} // namespace mxx

#endif // MXX_SAMPLESORT_HPP
