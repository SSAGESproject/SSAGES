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
 * @file    mpi_utils.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements some helpful MPI utility functions.
 */
#ifndef MXX_UTILS_HPP
#define MXX_UTILS_HPP

// MPI include
#include <mpi.h>

// C includes
#include <unistd.h> // for sleep()
#include <assert.h>

// C++ includes
#include <vector>
#include <algorithm>
#include <string>

// MXX
#include "datatypes.hpp"
#include "collective.hpp"


namespace mxx {

std::string get_processor_name() {
    char p_name[MPI_MAX_PROCESSOR_NAME+1];
    int p_len;
    MPI_Get_processor_name(p_name, &p_len);
    p_name[p_len] = '\0'; // make string NULL-terminated (if not yet so)
    std::string str(p_name);
    return str;
}

/**
 * @brief   Prints out summary information of which MPI processes (ranks)
 *          are running on which nodes.
 *
 * @param comm  The MPI communicator (default = MPI_COMM_WORLD)
 */
void print_node_distribution(MPI_Comm comm = MPI_COMM_WORLD)
{
    // get MPI parameters
    int rank;
    int p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);

    // get local processor name
    char p_name[MPI_MAX_PROCESSOR_NAME+1];
    int p_len;
    MPI_Get_processor_name(p_name, &p_len);
    p_name[p_len] = '\0'; // make string NULL-terminated (if not yet so)

    // gather all processor names to master
    // TODO: use gather with string specialization!
    std::vector<char> all_names_raw = gatherv(p_name,p_len+1, 0,comm);

    if (rank == 0)
    {
        std::vector<std::string> all_names(p);
        // disect the names into a set of strings
        std::vector<char>::iterator str_start = all_names_raw.begin();
        int i = 0;
        for (auto it = all_names_raw.begin(); it != all_names_raw.end(); ++it)
        {
            if (*it == '\0')
            {
                all_names[i++] = std::string(str_start, it);
                str_start = it+1;
            }
        }
        assert(i == p);
        std::map<std::string, std::vector<int> > procs_per_node;
        for (i = 0; i < p; ++i)
        {
            procs_per_node[all_names[i]].push_back(i);
        }
        // create array instead of map, then we can sort by first rank
        std::vector<std::pair<std::string, std::vector<int> > > proc_distr(procs_per_node.size());
        i = 0;
        for (auto it = procs_per_node.begin(); it != procs_per_node.end(); ++it)
        {
            // sort procs on each node
            std::sort(it->second.begin(), it->second.end());
            // put into the vector
            proc_distr[i++] = std::make_pair(it->first, it->second);
        }
        // sort the vector of node names by first rank
        std::sort(proc_distr.begin(), proc_distr.end(),
                  [](const std::pair<std::string, std::vector<int> >& x,
                     const std::pair<std::string, std::vector<int> >& y)
                  { return x.second.front() < y.second.front();});

        // print out the rank distribution
        std::cerr << "== Node distribution == " << std::endl;
        std::cerr << "== p=" << p << " processes on " << proc_distr.size() << " nodes ==" << std::endl;
        for (auto it = proc_distr.begin(); it != proc_distr.end(); ++it)
        {
            std::cerr << "--  Node: '" << it->first << "' (" << it->second.size() << "/" << p << ")" << std::endl;
            std::cerr << "        Ranks: ";
            for (auto rank_it = it->second.begin(); rank_it != it->second.end(); ++rank_it)
            {
                std::cerr << *rank_it;
                if (rank_it+1 != it->second.end())
                    std::cerr << ", ";
            }
            std::cerr << std::endl;
        }
    }
}

/**
 * @brief   Helper function for debugging. Prints out it's process PID
 *          and sleeps till the loop is broken via the debugger.
 *
 * @param wait_rank
 * @param comm
 */
void wait_gdb_attach(int wait_rank, MPI_Comm comm)
{
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    if (rank == wait_rank){
      std::cerr << "Rank " << rank << " is waiting in process " << getpid() << std::endl;
      int wait = 1;
      while (wait)
      {
        sleep(1);
      }
    }
    MPI_Barrier(comm);
}

} // namespace mxx

#endif // HPC_MPI_UTILS
