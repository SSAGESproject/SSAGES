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
 * @file    distribution.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements functions to redistribute data across processes.
 */

#ifndef MXX_DISTRIBUTION_HPP
#define MXX_DISTRIBUTION_HPP

#include <mpi.h>

#include <vector>
#include <queue>
#include "assert.h"

#include "datatypes.hpp"
#include "partition.hpp"
#include "reduction.hpp"
#include "collective.hpp"


#define MEASURE_LOAD_BALANCE 0
#if MEASURE_LOAD_BALANCE == 1
#include <cxx-prettyprint/prettyprint.hpp>
#endif

namespace mxx
{

/*
 * TODO:
 * - [x] new naming scheme (`(stable)_distribute_(inplace)` instead of `block_decomposition`)
 * - [x] iterator interfaces for implementation and most general versions
 * - [x] std::vector and std::basic_string specializations
 * - [ ] refactor the `partition.hpp` classes
 * - [x] get rid of code duplication
 * - [x] refactor API
 * - [ ] `redo_arbit_decomposition()` refactoring!
 * - [ ] surplus-send-pairing without allgather!
 */

namespace impl {
template<typename _InIterator, typename _OutIterator>
    void distribute_scatter(_InIterator begin, _InIterator end, _OutIterator out, size_t total_size, const mxx::comm& comm) {
        size_t local_size = std::distance(begin, end);
        // get root
        std::pair<size_t, int> max_pos = mxx::max_element(local_size, comm);
        int root = max_pos.second;
        MXX_ASSERT(mxx::all_same(root, comm));
        partition::block_decomposition<std::size_t> part(total_size, comm.size(), comm.rank());


        if (comm.rank() == root) {
            MXX_ASSERT(local_size == total_size);
            // use scatterv to send
            std::vector<size_t> send_sizes(comm.size());
            for (int i = 0; i < comm.size(); ++i) {
                send_sizes[i] = part.local_size(i);
            }
            // TODO: scatterv with iterators rather than pointers
            mxx::scatterv(&(*begin), send_sizes, &(*out), part.local_size(), root, comm);
        } else {
            // use scatterv to receive
            size_t recv_size = part.local_size();
            mxx::scatterv_recv(&(*out), recv_size, root, comm);
        }
    }
}

/**
 * @brief Fixes an unequal distribution into a block decomposition
 */
template<typename _InIterator, typename _OutIterator>
void stable_distribute(_InIterator begin, _InIterator end, _OutIterator out, const mxx::comm& comm) {
    // if there's only one process, return a copy
    if (comm.size() == 1) {
        std::copy(begin, end, out);
        return;
    }

    // get local and global size
    std::size_t local_size = std::distance(begin, end);
    std::size_t total_size = mxx::allreduce(local_size, comm);

    if (total_size == 0) {
        return;
    }

    // one process has all elements -> use scatter instead of all2all
    if (mxx::any_of(local_size == total_size, comm)) {
        impl::distribute_scatter(begin, end, out, total_size, comm);
    } else {
        // get prefix sum of size and total size
        std::size_t prefix = mxx::exscan(local_size, comm);

        // calculate where to send elements
        std::vector<size_t> send_counts(comm.size(), 0);
        partition::block_decomposition<std::size_t> part(total_size, comm.size(), comm.rank());
        int first_p = part.target_processor(prefix);
        std::size_t left_to_send = local_size;
        for (; left_to_send > 0 && first_p < comm.size(); ++first_p) {
            std::size_t nsend = std::min<std::size_t>(part.prefix_size(first_p) - prefix, left_to_send);
            send_counts[first_p] = nsend;
            left_to_send -= nsend;
            prefix += nsend;
        }
        std::vector<size_t> recv_counts = mxx::all2all(send_counts, comm);
        // TODO: accept iterators in mxx::all2all?
        mxx::all2allv(&(*begin), send_counts, &(*out), recv_counts, comm);
    }
}



namespace impl {
typedef std::make_signed<size_t>::type signed_size_t;

// negative `surpluses` represents a deficit
inline std::vector<size_t> surplus_send_pairing(std::vector<signed_size_t>& surpluses, int p, int rank, bool send_deficit = true) {
    // calculate the send and receive counts by a linear scan over
    // the surpluses, using a queue to keep track of all surpluses
    std::vector<size_t> send_counts(p, 0);
    std::queue<std::pair<int, signed_size_t> > fifo;
    for (int i = 0; i < p; ++i) {
        if (surpluses[i] == 0)
            continue;
        if (fifo.empty()) {
            fifo.push(std::make_pair(i, surpluses[i]));
        } else if (surpluses[i] > 0) {
            if (fifo.front().second > 0) {
                fifo.push(std::make_pair(i, surpluses[i]));
            } else {
                while (surpluses[i] > 0 && !fifo.empty())
                {
                    long long min = std::min(surpluses[i], -fifo.front().second);
                    int j = fifo.front().first;
                    surpluses[i] -= min;
                    fifo.front().second += min;
                    if (fifo.front().second == 0)
                        fifo.pop();
                    // these processors communicate!
                    if (rank == i)
                        send_counts[j] += min;
                    else if (rank == j && send_deficit)
                        send_counts[i] += min;
                }
                if (surpluses[i] > 0)
                    fifo.push(std::make_pair(i, surpluses[i]));
            }
        } else if (surpluses[i] < 0) {
            if (fifo.front().second < 0) {
                fifo.push(std::make_pair(i, surpluses[i]));
            } else {
                while (surpluses[i] < 0 && !fifo.empty())
                {
                    long long min = std::min(-surpluses[i], fifo.front().second);
                    int j = fifo.front().first;
                    surpluses[i] += min;
                    fifo.front().second -= min;
                    if (fifo.front().second == 0)
                        fifo.pop();
                    // these processors communicate!
                    if (rank == i && send_deficit)
                        send_counts[j] += min;
                    else if (rank == j)
                        send_counts[i] += min;
                }
                if (surpluses[i] < 0)
                    fifo.push(std::make_pair(i, surpluses[i]));
            }
        }
    }
    MXX_ASSERT(fifo.empty());

    return send_counts;
}
} // namespace impl


// non-stable version of `stable_distribute`
template<typename _InIterator, typename _OutIterator>
void distribute(_InIterator begin, _InIterator end, _OutIterator out, const mxx::comm& comm) {
    // if single process, copy input to output
    if (comm.size() == 1) {
        std::copy(begin, end, out);
        return;
    }

    // get sizes
    size_t local_size = std::distance(begin, end);
    size_t total_size = mxx::allreduce(local_size, comm);
    typedef typename std::iterator_traits<_InIterator>::value_type T;

    if (mxx::any_of(local_size == total_size, comm)) {
        // use scatterv instead of all2all based communication
        impl::distribute_scatter(begin, end, out, total_size, comm);
    } else {
        // use surplus send-pairing to minimize total communication volume
        // TODO: use all2all or send/recv depending on the maximum number of
        //       paired processes

        // get surplus
        partition::block_decomposition<std::size_t> part(total_size, comm.size(), comm.rank());
        impl::signed_size_t surplus = (impl::signed_size_t)local_size - (impl::signed_size_t)part.local_size();

        // allgather surpluses
        // TODO figure out how to do surplus send pairing without requiring allgather
        std::vector<impl::signed_size_t> surpluses = mxx::allgather(surplus, comm);
        MXX_ASSERT(std::accumulate(surpluses.begin(), surpluses.end(), static_cast<impl::signed_size_t>(0)) == 0);

        // get send counts
        std::vector<size_t> send_counts = impl::surplus_send_pairing(surpluses, comm.size(), comm.rank(), false);
        std::vector<size_t> recv_counts = mxx::all2all(send_counts, comm);

        if (surplus > 0) {
            MXX_ASSERT(std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0)) == 0);
            // TODO: use iterators not pointers
            mxx::all2allv(&(*(begin+((impl::signed_size_t)local_size-surplus))), send_counts, (T*)nullptr, recv_counts, comm);
            std::copy(begin, begin+(local_size-surplus), out);
        } else {
            MXX_ASSERT(std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0)) == 0);
            mxx::all2allv((const T*)nullptr, send_counts, &(*(out+local_size)), recv_counts, comm);
            std::copy(begin, end, out);
        }
    }
}

namespace impl {
template <class Container>
struct distribute_container {
    inline static Container stable_distribute(const Container& c, const mxx::comm& comm) {
        if (comm.size() == 1)
            return c;
        Container result;
        size_t local_size = c.size();
        // TODO: this allreduce is duplicated in the `stable_distribute` function
        //       which could be avoided
        size_t total_size = mxx::allreduce(local_size, comm);
        partition::block_decomposition<std::size_t> part(total_size, comm.size(), comm.rank());
        result.resize(part.local_size());

        // call the iterator based implementation
        ::mxx::stable_distribute(c.begin(), c.end(), result.begin(), comm);

        // return the container
        return result;
    }

    inline static void stable_distribute_inplace(Container& c, const mxx::comm& comm) {
        // we actually need a buffer anyway, so we use the non-inplace version
        // which returns a new container, and then assign it to the given container
        c = distribute_container::stable_distribute(c, comm);
    }

    inline static Container distribute(const Container& c, const mxx::comm& comm) {
        if (comm.size() == 1)
            return c;
        Container result;
        size_t local_size = c.size();
        // TODO: this allreduce is duplicated in the `distribute` function
        //       which could be avoided
        size_t total_size = mxx::allreduce(local_size, comm);
        partition::block_decomposition<std::size_t> part(total_size, comm.size(), comm.rank());
        result.resize(part.local_size());

        // call the iterator based implementation
        ::mxx::distribute(std::begin(c), std::end(c), std::begin(result), comm);

        // return the container
        return result;
    }

    inline static void distribute_inplace(Container& c, const mxx::comm& comm) {
        // custom implementation: don't need full sized buffers, only for receiving surpluses)
        // if single process, don't do anything
        if (comm.size() == 1)
            return;

        typedef typename std::iterator_traits<decltype(std::begin(c))>::value_type T;

        // get sizes
        size_t local_size = c.size();
        size_t total_size = mxx::allreduce(local_size, comm);
        partition::block_decomposition<std::size_t> part(total_size, comm.size(), comm.rank());

        if (mxx::any_of(local_size == total_size, comm)) {
            // use scatterv instead of all2all based communication
            Container result;
            result.resize(part.local_size());
            impl::distribute_scatter(std::begin(c), std::end(c), std::begin(result), total_size, comm);
            c.swap(result);
        } else {
            // use surplus send-pairing to minimize total communication volume
            // TODO: use all2all or send/recv depending on the maximum number of
            //       paired processes

            // get surplus
            impl::signed_size_t surplus = (impl::signed_size_t)local_size - (impl::signed_size_t)part.local_size();

            // allgather surpluses
            // TODO figure out how to do surplus send pairing without requiring allgather
            std::vector<impl::signed_size_t> surpluses = mxx::allgather(surplus, comm);
            MXX_ASSERT(std::accumulate(surpluses.begin(), surpluses.end(), static_cast<impl::signed_size_t>(0)) == 0);

            // get send counts
            std::vector<size_t> send_counts = impl::surplus_send_pairing(surpluses, comm.size(), comm.rank(), false);
            std::vector<size_t> recv_counts = mxx::all2all(send_counts, comm);

            if (surplus > 0) {
                MXX_ASSERT(std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0)) == 0);
                // TODO: use iterators not pointers
                mxx::all2allv(&(*(std::begin(c)+((impl::signed_size_t)local_size-surplus))), send_counts, (T*)nullptr, recv_counts, comm);
                c.resize(part.local_size());
            } else {
                MXX_ASSERT(std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0)) == 0);
                c.resize(part.local_size());
                mxx::all2allv((T*)nullptr, send_counts, &(*(std::begin(c)+local_size)), recv_counts, comm);
            }
        }
    }
};
} // namespace impl

/*
 * Re-distirbuted the vector into a perfect block decomposition.
 * (This invalidates all previous iterators)
 */
// std::vector specialization
template <typename T>
std::vector<T> stable_distribute(const std::vector<T>& local_els, const mxx::comm& comm) {
    return impl::distribute_container<std::vector<T>>::stable_distribute(local_els, comm);
}

// std::basic_string specialization
template <typename CharT, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT>>
std::basic_string<CharT, Traits, Alloc> stable_distribute(const std::basic_string<CharT, Traits, Alloc>& local_string, const mxx::comm& comm) {
    return impl::distribute_container<std::basic_string<CharT, Traits, Alloc>>::stable_distribute(local_string, comm);
}

// std::vector specialization
template <typename T>
void stable_distribute_inplace(std::vector<T>& local_els, const mxx::comm& comm) {
    impl::distribute_container<std::vector<T>>::stable_distribute_inplace(local_els, comm);
}

// std::basic_string specialization
template <typename CharT, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT>>
void stable_distribute_inplace(std::basic_string<CharT, Traits, Alloc>& local_string, const mxx::comm& comm) {
    impl::distribute_container<std::basic_string<CharT, Traits, Alloc>>::stable_distribute_inplace(local_string, comm);
}

// stable inplace

// non-stable variants:

// std::vector specialization
template <typename T>
std::vector<T> distribute(const std::vector<T>& local_els, const mxx::comm& comm) {
    return impl::distribute_container<std::vector<T>>::distribute(local_els, comm);
}

// std::basic_string specialization
template <typename CharT, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT>>
std::basic_string<CharT, Traits, Alloc> distribute(const std::basic_string<CharT, Traits, Alloc>& local_string, const mxx::comm& comm) {
    return impl::distribute_container<std::basic_string<CharT, Traits, Alloc>>::distribute(local_string, comm);
}

// non-stable inplace

// std::vector specialization
template <typename T>
void distribute_inplace(std::vector<T>& local_els, const mxx::comm& comm) {
    impl::distribute_container<std::vector<T>>::distribute_inplace(local_els, comm);
}

// std::basic_string specialization
template <typename CharT, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT>>
void distribute_inplace(std::basic_string<CharT, Traits, Alloc>& local_string, const mxx::comm& comm) {
    impl::distribute_container<std::basic_string<CharT, Traits, Alloc>>::distribute_inplace(local_string, comm);
}


/**
 * @brief Redistributes elements from the given decomposition across processors
 *        into the decomposition given by the requested local_size
 */
template<typename _InIterator, typename _OutIterator>
void redo_arbit_decomposition(_InIterator begin, _InIterator end, _OutIterator out, std::size_t new_local_size, const mxx::comm& comm) {
    // get local size
    std::size_t local_size = std::distance(begin, end);

    // if single process, simply copy to output
    if (comm.size() == 1) {
        MXX_ASSERT(local_size == new_local_size);
        std::copy(begin, end, out);
        return;
    }

    // get prefix sum of size and total size
#if !defined(NDEBUG) || MEASURE_LOAD_BALANCE != 0
    size_t total_size = mxx::allreduce(local_size, comm);
#endif
    size_t prefix = mxx::exscan(local_size, comm);
    if (comm.rank() == 0)
        prefix = 0;

#if MEASURE_LOAD_BALANCE
    size_t min = mxx::reduce(local_size, 0, mxx::min<size_t>(), comm);
    size_t max = mxx::reduce(local_size, 0, mxx::max<size_t>(), comm);
    size_t min_new = mxx::reduce(new_local_size, 0, mxx::min<size_t>(), comm);
    size_t max_new = mxx::reduce(new_local_size, 0, mxx::max<size_t>(), comm);
    if(comm.rank() == 0)
      std::cerr << " Decomposition: old [" << min << "," << max << "], new= [" << min_new << "," << max_new << "], for n=" << total_size << " fair decomposition: " << total_size / comm.size() << std::endl;

    std::vector<std::size_t> toReceive = mxx::gather(new_local_size, 0, comm);
    if(comm.rank() == 0)
      std::cerr << toReceive << std::endl;
#endif

    // get the new local sizes from all processors
    // NOTE: this all-gather is what makes the arbitrary decomposition worse
    // in terms of complexity than when assuming a block decomposition
    std::vector<std::size_t> new_local_sizes = mxx::allgather(new_local_size, comm);
    MXX_ASSERT(std::accumulate(new_local_sizes.begin(), new_local_sizes.end(), static_cast<size_t>(0)) == total_size);

    // calculate where to send elements
    std::vector<size_t> send_counts(comm.size(), 0);
    int first_p;
    std::size_t new_prefix = 0;
    for (first_p = 0; first_p < comm.size()-1; ++first_p)
    {
        // find processor for which the prefix sum exceeds mine
        // i have to send to the previous
        if (new_prefix + new_local_sizes[first_p] > prefix)
            break;
        new_prefix += new_local_sizes[first_p];
    }

    //= block_partition_target_processor(total_size, p, prefix);
    std::size_t left_to_send = local_size;
    for (; left_to_send > 0 && first_p < comm.size(); ++first_p)
    {
        // make the `new` prefix inclusive (is an exlcusive prefix prior)
        new_prefix += new_local_sizes[first_p];
        // send as many elements to the current processor as it needs to fill
        // up, but at most as many as I have left
        std::size_t nsend = std::min<std::size_t>(new_prefix - prefix, left_to_send);
        send_counts[first_p] = nsend;
        // update the number of elements i have left (`left_to_send`) and
        // at which global index they start `prefix`
        left_to_send -= nsend;
        prefix += nsend;
    }

    std::vector<size_t> recv_counts = mxx::all2all(send_counts, comm);
    // TODO: all2allv for iterators
    mxx::all2allv(&(*begin), send_counts, &(*out), recv_counts, comm);
}

// function for block decompsing vector of two partitions (for equal
// distributing of two halves)
// this one is stable for both halves.
template<typename _Iterator>
_Iterator stable_block_decompose_partitions(_Iterator begin, _Iterator mid, _Iterator end, const mxx::comm& comm)
{
    typedef typename std::iterator_traits<_Iterator>::value_type T;
    // return same sequence if there's a single process
    if (comm.size() == 1)
        return mid;

    // get sizes
    std::size_t left_local_size = std::distance(begin, mid);
    std::size_t right_local_size = std::distance(mid, end);
    // TODO: use array of 2 (single reduction!)
    std::size_t left_size = mxx::allreduce(left_local_size, comm);
    std::size_t right_size = mxx::allreduce(right_local_size, comm);
    partition::block_decomposition<std::size_t> part(left_size+right_size, comm.size(), comm.rank());
    partition::block_decomposition<std::size_t> left_part(left_size, comm.size(), comm.rank());

    // shuffle into buffer
    std::vector<T> buffer(part.local_size());
    redo_block_decomposition(begin, mid, buffer.begin(), comm);
    redo_arbit_decomposition(mid, end, buffer.begin()+left_part.local_size(), part.local_size() - left_part.local_size(), comm);

    // copy back
    std::copy(buffer.begin(), buffer.end(), begin);
    return begin + left_part.local_size();
}

// non-stable version (much faster, since data exchange is only for unequal parts)
template<typename _Iterator>
_Iterator block_decompose_partitions(_Iterator begin, _Iterator mid, _Iterator end, const mxx::comm& comm)
{
    typedef typename std::iterator_traits<_Iterator>::value_type T;
    // return same sequence if there's a single process
    if (comm.size() == 1)
        return mid;
    // get sizes
    size_t left_local_size = std::distance(begin, mid);
    size_t left_size = mxx::allreduce(left_local_size, comm);
    partition::block_decomposition<std::size_t> left_part(left_size, comm.size(), comm.rank());
    impl::signed_size_t surplus = (impl::signed_size_t)left_local_size - (impl::signed_size_t)left_part.local_size();
    bool fits = end-begin >= left_part.local_size();
    bool all_fits = mxx::all_of(fits, comm);
    if (!all_fits)
        return mid;
    std::vector<impl::signed_size_t> surpluses = mxx::allgather(surplus, comm);

    MXX_ASSERT(std::accumulate(surpluses.begin(), surpluses.end(), static_cast<size_t>(0)) == 0);

    // get send counts
    std::vector<size_t> send_counts = impl::surplus_send_pairing(surpluses, comm.size(), comm.rank(), true);

    std::size_t send_size = std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0));
    std::vector<T> buffer;
    if (send_size > 0)
        buffer.resize(send_size);

    std::vector<size_t> recv_counts = mxx::all2all(send_counts, comm);
    // assert send and receive size are the same!
    for (int i = 0; i < comm.size(); ++i) {
        MXX_ASSERT(send_counts[i] == recv_counts[i]);
    }

    // send from the surplus, receive into buffer
    if (surplus > 0) {
        mxx::all2allv(&(*(mid - surplus)), send_counts, &buffer[0], recv_counts, comm);
        std::copy(buffer.begin(), buffer.end(), mid-surplus);
    } else if (surplus < 0) {
        mxx::all2allv(&(*mid), send_counts, &buffer[0], recv_counts, comm);
        std::copy(buffer.begin(), buffer.end(), mid);
    } else {
        MXX_ASSERT(send_size == 0);
        mxx::all2allv((T*)nullptr, send_counts, (T*)nullptr, recv_counts, comm);
    }

    // copy back
    return mid - surplus;
}
} // namespace mxx

#endif // MXX_DISTRIBUTION_HPP
