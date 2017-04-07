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
 * @file    bcast.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @group   collective
 * @brief   Bcast operations.
 *
 */

#ifndef MXX_BCAST_HPP
#define MXX_BCAST_HPP

#include <mpi.h>
#include <vector>
#include <string>

// mxx includes
#include "common.hpp"
#include "datatypes.hpp"
#include "comm_fwd.hpp"
#include "reduction.hpp"

namespace mxx {

/// Generic bcast
template <typename T>
void bcast(T* data, size_t count, int root, const mxx::comm& comm) {
    mxx::datatype dt = mxx::get_datatype<T>();
    if (count >= mxx::max_int) {
        mxx::datatype bigdt = dt.contiguous(count);
        MPI_Bcast(data, 1, bigdt.type(), root, comm);
    } else {
        MPI_Bcast(data, count, dt.type(), root, comm);
    }
}

/// Broadcast single value
template <typename T>
void bcast(T& value, int root, const mxx::comm& comm) {
    mxx::datatype dt = mxx::get_datatype<T>();
    MPI_Bcast(&value, 1, dt.type(), root, comm);
}

// container bcast implementation
namespace impl {
template <typename Container>
struct container_bcast {
static void do_bcast(Container& c, int root, const mxx::comm& comm) {
    // first bcast the size, so that receiving can allocate the string
    size_t size = c.size();
    mxx::bcast(size, root, comm);
    if (comm.rank() != root)
        c.resize(size);
    mxx::bcast(&c[0], size, root, comm);
}
};

} // namespace impl

// broadcast a vector
template <typename T, typename Alloc = std::allocator<T>>
void bcast(std::vector<T, Alloc>& vec, int root, const mxx::comm& comm) {
    impl::container_bcast<std::vector<T, Alloc>>::do_bcast(vec, root, comm);
}

// broadcast a string
template <typename CharT = char, typename Traits = std::char_traits<CharT>,
          typename Alloc = std::allocator<CharT>>
void bcast(std::basic_string<CharT, Traits, Alloc>& str, int root, const mxx::comm& comm) {
    impl::container_bcast<std::basic_string<CharT, Traits, Alloc>>::do_bcast(str, root, comm);
}

} // namespace mxx

#endif // MXX_BCAST_HPP
