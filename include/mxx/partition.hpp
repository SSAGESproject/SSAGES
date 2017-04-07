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
 * @file    paritition.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the partition API and different data partitions, most
 *          notably the block decomposition.
 */

#ifndef MXX_PARTITION_HPP
#define MXX_PARTITION_HPP

#include <vector>
#include "assert.h"

namespace mxx
{
namespace partition
{

// We want inlined non virtual functions for the partition. We thus need to use
// templating when these classes are used and no real inheritance to enforce
// the API.
template <typename index_t>
class block_decomposition
{
public:
    block_decomposition() : n(0), p(0), rank(0) {}

    /// Constructor (no default construction)
    block_decomposition(index_t n, int p, int rank)
        : n(n), p(p), rank(rank)
    {
    }

    index_t local_size()
    {
        return n/p + ((static_cast<index_t>(rank) < (n % static_cast<index_t>(p))) ? 1 : 0);
    }

    index_t local_size(int rank)
    {
        return n/p + ((static_cast<index_t>(rank) < (n % static_cast<index_t>(p))) ? 1 : 0);
    }

    index_t prefix_size()
    {
        return (n/p)*(rank+1) + std::min<index_t>(n % p, rank + 1);
    }

    index_t prefix_size(int rank)
    {
        return (n/p)*(rank+1) + std::min<index_t>(n % p, rank + 1);
    }

    index_t excl_prefix_size()
    {
        return (n/p)*rank + std::min<index_t>(n % p, rank);
    }

    index_t excl_prefix_size(int rank)
    {
        return (n/p)*rank + std::min<index_t>(n % p, rank);
    }

    // which processor the element with the given global index belongs to
    int target_processor(index_t global_index)
    {
        if (global_index < ((n/p)+1)*(n % p))
        {
            // a_i is within the first n % p processors
            return global_index/((n/p)+1);
        }
        else
        {
            return n%p + (global_index - ((n/p)+1)*(n % p))/(n/p);
        }
    }

    /// Destructor
    virtual ~block_decomposition () {}
private:
    /* data */
    /// Number of elements
    index_t n;
    /// Number of processors
    int p;
    /// Processor rank
    int rank;
};


template <typename index_t>
class block_decomposition_buffered
{
public:
    block_decomposition_buffered() {}

    block_decomposition_buffered(index_t n, int p, int rank)
        : n(n), p(p), rank(rank), div(n / p), mod(n % p),
          loc_size(div + (static_cast<index_t>(rank) < mod ? 1 : 0)),
          prefix(div*rank + std::min<index_t>(mod, rank)),
          div1mod((div+1)*mod)
    {
    }

    block_decomposition_buffered(const block_decomposition_buffered& other)
        : n(other.n), p(other.p), rank(other.rank), div(other.div),
        mod(other.mod), loc_size(other.loc_size), prefix(other.prefix),
        div1mod(other.div1mod) {}

    block_decomposition_buffered& operator=(const block_decomposition_buffered& other) = default;

    index_t local_size()
    {
        return loc_size;
    }

    index_t local_size(int rank)
    {
        return div + (static_cast<index_t>(rank) < mod ? 1 : 0);
    }

    index_t prefix_size()
    {
        return prefix + loc_size;
    }

    index_t prefix_size(int rank)
    {
        return div*(rank+1) + std::min<index_t>(mod, rank + 1);
    }

    index_t excl_prefix_size()
    {
        return prefix;
    }

    index_t excl_prefix_size(int rank)
    {
        return div*rank + std::min<index_t>(mod, rank);
    }

    // which processor the element with the given global index belongs to
    int target_processor(index_t global_index)
    {
        // TODO: maybe also buffer (div+1)*mod, would save one multiplication
        // in each call to this
        if (global_index < div1mod)
        {
            // a_i is within the first n % p processors
            return global_index/(div+1);
        }
        else
        {
            return mod + (global_index - div1mod)/div;
        }
    }

    virtual ~block_decomposition_buffered () {}
private:
    /* data */
    /// Number of elements
    index_t n;
    /// Number of processors
    int p;
    /// Processor rank
    int rank;

    // derived/buffered values (for faster computation of results)
    index_t div; // = n/p
    index_t mod; // = n%p
    // local size (number of local elements)
    index_t loc_size;
    // the exclusive prefix (number of elements on previous processors)
    index_t prefix;
    /// number of elements on processors with one more element
    index_t div1mod; // = (n/p + 1)*(n % p)
};
} // namespace partition

/**
 * @brief Returns a block partitioning of an input of size `n` among `p` processors.
 *
 * @param n The number of elements.
 * @param p The number of processors.
 *
 * @return A vector of the number of elements for each processor.
 */
/*
std::vector<int> block_partition(int n, int p)
{
    // init result
    std::vector<int> partition(p);
    // get the number of elements per processor
    int local_size = n / p;
    // and the elements that are not evenly distributable
    int remaining = n % p;
    for (int i = 0; i < p; ++i) {
        if (i < remaining) {
            partition[i] = local_size + 1;
        } else {
            partition[i] = local_size;
        }
    }
    return partition;
}
*/

} // namespace mxx

#endif // MXX_PARTITION_HPP
