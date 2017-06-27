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
 * @file    timer.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   MPI synchronized section timer.
 */

#ifndef MXX_TIMER_HPP
#define MXX_TIMER_HPP

// MPI include
#include <mpi.h>

// C++ includes
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>

// mxx includes
#include "reduction.hpp"

namespace mxx
{

template <typename duration>
class timer_impl {
protected:
    // use the monotonic `steady_clock` for time measurements of code
    std::chrono::steady_clock::time_point start;
public:
    // constructor
    timer_impl() {
        start = std::chrono::steady_clock::now();
    }

    // returns time elapsed since creation
    typename duration::rep elapsed() const {
      std::chrono::steady_clock::time_point stop = std::chrono::steady_clock::now();
      typename duration::rep elapsed_time = duration(stop-start).count();
      return elapsed_time;
    }
};

/// mxx::timer: Specialization for measuring milliseconds in double precision
using timer = timer_impl<std::chrono::duration<double, std::milli> >;

template <typename duration>
class section_timer_impl {
protected:
    /// The output stream
    std::ostream& ostr;
    /// The MPI communicator for barrier sync and reductions
    MPI_Comm comm;
    /// the root process (only the process with this rank prints the timing
    /// results)
    int root;
    /// The seperator used for the machine readable output (default: '\t')
    std::string sep;
    /// The current depth of timer instantiations.
    static int depth;
    /// The type of the time duration
    typedef typename duration::rep time_rep;

    // use the monotonic `steady_clock` for time measurements of code
    /// Type of a time point
    std::chrono::steady_clock::time_point start;

public:
    /// Construtor with reasonable defaults
    section_timer_impl(std::ostream& outstream = std::cerr, MPI_Comm comm = MPI_COMM_WORLD, int root = 0)
      : ostr(outstream), comm(comm), root(root), sep("\t") {
        // sync on barrier
        MPI_Barrier(comm);
        start = std::chrono::steady_clock::now();
        section_timer_impl::depth++;
    }

    /// End the section with the given name. This prints out the elapsed time
    /// information aggregated from all processors
    void end_section(const std::string& name) {
        std::chrono::steady_clock::time_point stop = std::chrono::steady_clock::now();
        time_rep elapsed = duration(stop-start).count();
        // TODO: reduce rather than allreduce
        // TODO: use single reduction with custom operator
        double delapsed = elapsed;
        double min = mxx::allreduce(delapsed, mxx::min<double>(), comm);
        double max = mxx::allreduce(delapsed, mxx::max<double>(), comm);
        double sum = mxx::allreduce(delapsed, std::plus<time_rep>(), comm);
        // calc mean:
        int p, rank;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &p);
        double mean = sum / (double)p;
        // only root process outputs the timings
        if (rank == root)
            ostr << std::setprecision(3) << std::scientific << "TIMER" << sep << min << sep << mean << sep << max << sep << depth << sep << name << std::endl;
        // restart timer
        MPI_Barrier(comm);
        start = std::chrono::steady_clock::now();
    }

    virtual ~section_timer_impl() {
      section_timer_impl::depth--;
    }
};

class empty_section_timer_impl {
    public:
    empty_section_timer_impl(std::ostream& = std::cerr, MPI_Comm = MPI_COMM_WORLD, int = 0) {
    }

    inline void end_section(const std::string&) {
    }
};

template <typename duration>
int section_timer_impl<duration>::depth = 0;

/// mxx::section_timer: specialization for measuring milliseconds in double precision
#ifdef MXX_DISABLE_TIMER
using section_timer = empty_section_timer_impl;
#else
using section_timer = section_timer_impl<std::chrono::duration<double, std::milli> >;
#endif

} // namespace mxx



#endif // MXX_TIMER_HPP
