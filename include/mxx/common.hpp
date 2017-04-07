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

#ifndef MXX_COMMON_HPP
#define MXX_COMMON_HPP

#include <mpi.h>
#include <limits>

namespace mxx {

// define a max limit (either as MAX_INT or as compilation defined variable)
#ifdef MXX_MAX_INT
// set to smaller value for testing
constexpr size_t max_int = MXX_MAX_INT;
#else
/// maximum message size for MPI
constexpr size_t max_int = std::numeric_limits<int>::max();
#endif

// check type of MPI_Aint
static_assert(sizeof(size_t) == sizeof(MPI_Aint), "MPI_Aint must be the same size as size_t");

// assertion failure
inline void assert_fail(const char * cond, const char* file, int line, const char* func) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("\n[Rank %d] %s:%d: %s: Assertion `%s` failed. Aborting.\n", rank, file, line, func, cond);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, -1);
}

} // namespace mxx


// assert macro
#ifdef NDEBUG
#define MXX_ASSERT(cond) ((void)0) /* TODO should we really deactivate all asserts in a -DNDEBUG build? */
#elif defined(MXX_NDEBUG)
#define MXX_ASSERT(cond) ((void)0)
#else
#define MXX_ASSERT(cond) {if (!(cond)) { mxx::assert_fail("" #cond "",__FILE__, __LINE__, __func__); }}
#endif

#endif // MXX_COMMON_HPP
