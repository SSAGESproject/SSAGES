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
 * @file    samplesort.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements parallel, MPI sample sort
 */

#ifndef MXX_SAMPLESORT_HPP
#define MXX_SAMPLESORT_HPP

#include "comm.hpp"
#include "partition.hpp"
#include "datatypes.hpp"
#include "collective.hpp"
#include "shift.hpp"
#include "distribution.hpp"
#include "bitonicsort.hpp"

#include <iterator>
#include <algorithm>
#include <vector>

#include <cxx-prettyprint/prettyprint.hpp>

#ifdef __GNUC__
#ifndef __clang__
// for multiway-merge
// TODO: implement multiway merging ourselves
#include <parallel/multiway_merge.h>
#include <parallel/merge.h>
#define MXX_USE_GCC_MULTIWAY_MERGE
#endif
#endif


#define SS_ENABLE_TIMER 0
#if SS_ENABLE_TIMER
#include "timer.hpp"
#define SS_TIMER_START(comm) mxx::section_timer timer(std::cerr, comm, 0);
#define SS_TIMER_END_SECTION(str) timer.end_section(str);
#else
#define SS_TIMER_START(comm)
#define SS_TIMER_END_SECTION(str)
#endif

namespace mxx {
namespace impl {

template<typename _Iterator, typename _Compare>
bool is_sorted(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm)
{
    // get value type of underlying data
    typedef typename std::iterator_traits<_Iterator>::value_type value_type;

    // check that it is locally sorted (int for MPI_Reduction)
    bool sorted = std::is_sorted(begin, end, comp);

    // if single process, this is it
    if (comm.size() == 1)
        return sorted;

    // compare if last element on left processor is not bigger than first
    // element on mine
    value_type left_el = mxx::right_shift(*(end-1), comm);

    // check if sorted with respect to neighbors
    if (comm.rank() > 0) {
        sorted = sorted && !comp(*begin, left_el);
    }

    // now check that this condition is true for all processes
    return mxx::all_of(sorted, comm);
}

template <typename _Iterator, typename _Compare>
std::vector<typename std::iterator_traits<_Iterator>::value_type>
sample_arbit_decomp(_Iterator begin, _Iterator end, _Compare comp, size_t s, const mxx::comm& comm) {
    typedef typename std::iterator_traits<_Iterator>::value_type value_type;
    std::size_t local_size = std::distance(begin, end);

    // get total size n
    std::size_t total_size = mxx::allreduce(local_size, comm);

    // p = number of processes
    int p = comm.size();

    //  pick a total of s*p samples, thus locally pick ceil((local_size/n)*s*p)
    //  and at least one samples from each processor.
    //  this will result in at least s*p samples.
    std::size_t local_s;
    if (local_size == 0)
        local_s = 0;
    else
        local_s = std::max<std::size_t>(((local_size*s*p)+total_size-1)/total_size, 1);

    MXX_ASSERT(mxx::allreduce(local_s, comm) >= p*s);

    // init samples
    std::vector<value_type> local_splitters;

    // pick local samples
    if (local_s > 0) {
        local_splitters.resize(local_s);
        _Iterator pos = begin;
        for (std::size_t i = 0; i < local_splitters.size(); ++i)
        {
            std::size_t bucket_size = local_size / (local_s+1) + (i < (local_size % (local_s+1)) ? 1 : 0);
            // pick last element of each bucket
            pos += (bucket_size-1);
            local_splitters[i] = *pos;
            ++pos;
        }
    }

    // distribute elements equally
    mxx::distribute_inplace(local_splitters, comm);

    // discard extra splitters, to make it even
    if (local_splitters.size() != s) {
        MXX_ASSERT(local_splitters.size() > s);
        local_splitters.resize(s);
    }

    // sort splitters using parallel bitonic sort
    bitonic_sort(local_splitters.begin(), local_splitters.end(), comp, comm);

    // select the last element on each process but the last
    value_type my_splitter = local_splitters.back();
    // allgather splitters (from all but the last processor)
    std::vector<size_t> recv_sizes(comm.size(), 1);
    recv_sizes.back() = 0;
    std::vector<value_type> result_splitters = mxx::allgatherv(&my_splitter, (comm.rank() + 1 < comm.size()) ? 1 : 0, recv_sizes, comm);

    // return resulting splitters
    return result_splitters;
}


template <typename _Iterator, typename _Compare>
std::vector<typename std::iterator_traits<_Iterator>::value_type>
sample_block_decomp(_Iterator begin, _Iterator end, _Compare comp, size_t s, const mxx::comm& comm) {
    typedef typename std::iterator_traits<_Iterator>::value_type value_type;
    std::size_t local_size = std::distance(begin, end);
    MXX_ASSERT(local_size > 0);

    // 1. samples
    //  - pick `s` samples equally spaced such that `s` samples define `s+1`
    //    subsequences in the sorted order
    std::vector<value_type> local_splitters(s);
    _Iterator pos = begin;
    for (std::size_t i = 0; i < local_splitters.size(); ++i) {
        std::size_t bucket_size = local_size / (s+1) + (i < (local_size % (s+1)) ? 1 : 0);
        // pick last element of each bucket
        pos += (bucket_size-1);
        local_splitters[i] = *pos;
        ++pos;
    }

    // sort splitters using parallel bitonic sort
    bitonic_sort(local_splitters.begin(), local_splitters.end(), comp, comm);

    // select the last element on each process but the last
    value_type my_splitter = local_splitters.back();
    // allgather splitters (from all but the last processor)
    std::vector<size_t> recv_sizes(comm.size(), 1);
    recv_sizes.back() = 0;
    std::vector<value_type> result_splitters = mxx::allgatherv(&my_splitter, (comm.rank() + 1 < comm.size()) ? 1 : 0, recv_sizes, comm);

    // return resulting splitters
    return result_splitters;
}

template <typename _Iterator, typename _Compare>
std::vector<size_t> split(_Iterator begin, _Iterator end, _Compare comp, const std::vector<typename std::iterator_traits<_Iterator>::value_type>& splitters, const mxx::comm& comm) {
    // 5. locally find splitter positions in data
    //    (if an identical splitter appears at least three times (or more),
    //    then split the intermediary buckets evenly) => send_counts
    MXX_ASSERT(splitters.size() == (size_t)comm.size() - 1);
    std::vector<size_t> send_counts(comm.size(), 0);
    _Iterator pos = begin;
    size_t local_size = std::distance(begin, end);
    partition::block_decomposition<std::size_t> local_part(local_size, comm.size(), comm.rank());
    for (std::size_t i = 0; i < splitters.size();) {
        // get the number of splitters which are equal starting from `i`
        unsigned int split_by = 1;
        while (i+split_by < splitters.size()
               && !comp(splitters[i], splitters[i+split_by])) {
            ++split_by;
        }

        // get the range of equal elements
        std::pair<_Iterator, _Iterator> eqr = std::equal_range(pos, end, splitters[i], comp);

        // assign smaller elements to processor left of splitter (= `i`)
        send_counts[i] += std::distance(pos, eqr.first);
        pos = eqr.first;

        // split equal elements fairly across processors
        std::size_t eq_size = std::distance(pos, eqr.second);
        // try to split approx equal:
        std::size_t eq_size_split = (eq_size + send_counts[i]) / (split_by+1) + 1;
        for (unsigned int j = 0; j < split_by; ++j) {
            std::size_t out_size = 0;
            if (send_counts[i+j] < local_part.local_size(i+j)) {
                // try to distribute fairly
                out_size = std::min(std::max(local_part.local_size(i+j) - send_counts[i+j], eq_size_split), eq_size);
                eq_size -= out_size;
            }
            send_counts[i+j] += out_size;
        }
        // assign remaining elements to next processor
        send_counts[i+split_by] += eq_size;
        i += split_by;
        pos = eqr.second;
    }
    // send last elements to last processor
    std::size_t out_size = std::distance(pos, end);
    send_counts[comm.size() - 1] += out_size;
    MXX_ASSERT(std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0)) == local_size);
    return send_counts;
}

template <typename _Iterator, typename _Compare>
std::vector<size_t> stable_split(_Iterator begin, _Iterator end, _Compare comp, const std::vector<typename std::iterator_traits<_Iterator>::value_type>& splitters, const mxx::comm& comm) {
    // 5. locally find splitter positions in data
    //    (if an identical splitter appears at least three times (or more),
    //    then split the intermediary buckets evenly) => send_counts
    MXX_ASSERT(splitters.size() == (size_t) comm.size() - 1);
    std::vector<size_t> send_counts(comm.size(), 0);
    _Iterator pos = begin;
    size_t local_size = std::distance(begin, end);
    partition::block_decomposition<std::size_t> local_part(local_size, comm.size(), comm.rank());
    for (std::size_t i = 0; i < splitters.size();) {
        // get the number of splitters which are equal starting from `i`
        unsigned int split_by = 1;
        while (i+split_by < splitters.size()
               && !comp(splitters[i], splitters[i+split_by])) {
            ++split_by;
        }

        // get the range of equal elements
        std::pair<_Iterator, _Iterator> eqr = std::equal_range(pos, end, splitters[i], comp);

        // assign smaller elements to processor left of splitter (= `i`)
        send_counts[i] += std::distance(pos, eqr.first);
        pos = eqr.first;

        // split equal elements fairly across processors
        std::size_t eq_size = std::distance(pos, eqr.second);

        // send equal elements to processor based on my own rank compared to
        // how many equal splitters there are
        if (split_by == 1) {
            // 1) if there is only one splitter, assign equal elements to next processor (no splitting)
            send_counts[i+1] += eq_size;
        } else {
            // 2) if there is >= 2 equal splitters: split processors into `split_by` regions
            unsigned int targetp = (comm.rank() * split_by) / comm.size();
            if (targetp >= split_by) targetp = split_by-1;
            send_counts[i+1+targetp] += eq_size;
        }
        i += split_by;
        pos = eqr.second;
    }

    // send last elements to last processor
    std::size_t out_size = std::distance(pos, end);
    send_counts[comm.size() - 1] += out_size;
    MXX_ASSERT(std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0)) == local_size);
    return send_counts;
}


template<typename _Iterator, typename _Compare, bool _Stable = false>
void samplesort(_Iterator begin, _Iterator end, _Compare comp, const mxx::comm& comm) {
    // get value type of underlying data
    typedef typename std::iterator_traits<_Iterator>::value_type value_type;

    int p = comm.size();

    SS_TIMER_START(comm);

    // perform local (stable) sorting
    if (_Stable)
        std::stable_sort(begin, end, comp);
    else
        std::sort(begin, end, comp);

    // sequential case: we're done
    if (p == 1)
        return;

#if SS_ENABLE_TIMER
    MPI_Barrier(comm);
#endif
    SS_TIMER_END_SECTION("local_sort");


    // local size
    std::size_t local_size = std::distance(begin, end);

    // check if we have a perfect block decomposition
    std::size_t global_size = mxx::allreduce(local_size, comm);
    partition::block_decomposition<std::size_t> mypart(global_size, p, comm.rank());
    bool _AssumeBlockDecomp = mxx::all_of(local_size == mypart.local_size(), comm);

    // sample sort
    // 1. local sort
    // 2. pick `s` samples regularly spaced on each processor
    // 3. bitonic sort samples
    // 4. allgather the last sample of each process -> splitters
    // 5. locally find splitter positions in data
    //    (if an identical splitter appears twice, then split evenly)
    //    => send_counts
    // 6. distribute send_counts with all2all to get recv_counts
    // 7. allocate enough space (may be more than previously allocated) for receiving
    // 8. all2allv
    // 9. local reordering (multiway-merge or again std::sort)
    // A. equalizing distribution into original size (e.g.,block decomposition)
    //    by sending elements to neighbors

    // get splitters, using the method depending on whether the input consists
    // of arbitrary decompositions or not
    std::vector<value_type> local_splitters;
    // number of samples
    size_t s = p-1;
    if(_AssumeBlockDecomp)
        local_splitters = sample_block_decomp(begin, end, comp, s, comm);
    else
        local_splitters = sample_arbit_decomp(begin, end, comp, s, comm);
    SS_TIMER_END_SECTION("get_splitters");

    // 5. locally find splitter positions in data
    //    (if an identical splitter appears at least three times (or more),
    //    then split the intermediary buckets evenly) => send_counts
    std::vector<size_t> send_counts;
    if (_Stable)
        send_counts = stable_split(begin, end, comp, local_splitters, comm);
    else
        send_counts = split(begin, end, comp, local_splitters, comm);
    SS_TIMER_END_SECTION("send_counts");


    std::vector<size_t> recv_counts = mxx::all2all(send_counts, comm);
    std::size_t recv_n = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));
    MXX_ASSERT(!_AssumeBlockDecomp || (local_size <= (size_t)p || recv_n <= 2* local_size));
    std::vector<value_type> recv_elements(recv_n);
    // TODO: use collective with iterators [begin,end) instead of pointers!
    mxx::all2allv(&(*begin), send_counts, &(*recv_elements.begin()), recv_counts, comm);
    SS_TIMER_END_SECTION("all2all");

    // 9. local reordering
    /* multiway-merge (using the implementation in __gnu_parallel) */
#ifdef MXX_USE_GCC_MULTIWAY_MERGE
    if (!_Stable && local_size > (size_t)p*p) {
        // prepare the sequence offsets
        typedef typename std::vector<value_type>::iterator val_it;
        std::vector<std::pair<val_it, val_it> > seqs(p);
        std::vector<size_t> recv_displs = mxx::impl::get_displacements(recv_counts);
        for (int i = 0; i < p; ++i) {
            seqs[i].first = recv_elements.begin() + recv_displs[i];
            seqs[i].second = seqs[i].first + recv_counts[i];
        }
        val_it start_merge_it = recv_elements.begin();

        std::size_t merge_n = local_size;
        value_type* merge_buf_begin = &(*begin);
        std::vector<value_type> merge_buf;
        // TODO: reasonable values for the buffer?
        // currently: at least 1/10 th of the size to merge or 1 MiB
        if (local_size == 0 || local_size < recv_n / 10)
        {
            // at least 1MB buffer
            merge_n = std::max<std::size_t>(recv_n / 10, (1024*1024)/sizeof(value_type));
            merge_buf.resize(merge_n);
            merge_buf_begin = &merge_buf[0];
        }
        for (; recv_n > 0;)
        {
            if (recv_n < merge_n)
                merge_n = recv_n;
            // i)   merge at most `local_size` many elements sequentially
            __gnu_parallel::sequential_tag seq_tag;
            __gnu_parallel::multiway_merge(seqs.begin(), seqs.end(), merge_buf_begin, merge_n, comp, seq_tag);

            // ii)  compact the remaining elements in `recv_elements`
            for (int i = p-1; i > 0; --i)
            {
                seqs[i-1].first = std::copy_backward(seqs[i-1].first, seqs[i-1].second, seqs[i].first);
                seqs[i-1].second = seqs[i].first;
            }
            // iii) copy the output buffer `local_size` elements back into
            //      `recv_elements`
            start_merge_it = std::copy(merge_buf_begin, merge_buf_begin + merge_n, start_merge_it);
            assert(start_merge_it == seqs[0].first);

            // reduce the number of elements to be merged
            recv_n -= merge_n;
        }
        // clean up
        merge_buf.clear(); merge_buf.shrink_to_fit();
    } else
#endif
    {
        if (_Stable)
            std::stable_sort(recv_elements.begin(), recv_elements.end(), comp);
        else
            std::sort(recv_elements.begin(), recv_elements.end(), comp);
    }

#if SS_ENABLE_TIMER
    MPI_Barrier(comm);
#endif
    SS_TIMER_END_SECTION("local_merge");

    // A. equalizing distribution into original size (e.g.,block decomposition)
    //    by elements to neighbors
    //    and save elements into the original iterator positions
    if (_AssumeBlockDecomp)
        stable_distribute(recv_elements.begin(), recv_elements.end(), begin, comm);
    else
        redo_arbit_decomposition(recv_elements.begin(), recv_elements.end(), begin, local_size, comm);

    SS_TIMER_END_SECTION("fix_partition");
}

} // namespace impl
} // namespace mxx

#endif // MXX_SAMPLESORT_HPP


