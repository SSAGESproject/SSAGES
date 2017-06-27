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

#ifndef MXX_COMM_FWD_HPP
#error "Never include this file directly"
#endif

// only include this file once all three dependencies are fully included
// i.e. the last of the three will include this file

#include "shift.hpp"
#ifdef MXX_SHIFT_DONE
#include "reduction.hpp"
#ifdef MXX_REDUCTION_DONE
#include "collective.hpp"
#ifdef MXX_COLLECTIVE_DONE
#include "stream.hpp"
#ifdef MXX_STREAM_DONE

#ifndef MXX_COMM_DEF_HPP
#define MXX_COMM_DEF_HPP

namespace mxx {

template <typename Func>
inline void comm::with_subset(bool cond, Func f) const {
    // only split communicator if the condition is not true everywhere
    // since splitting of communicators is an expensive operation
    if (!mxx::all_of(cond, *this)) {
        comm s(this->split(cond));
        if (cond) {
           f(s);
        }
    } else {
        f(*this);
    }
}

inline mxx::sync_ostream comm::sync_cout() const {
    return mxx::sync_cout(*this);
}

inline mxx::sync_ostream comm::sync_cerr() const {
    return mxx::sync_cerr(*this);
}

inline comm comm::split_shared() const {
        comm o;
#if MPI_VERSION >= 3
        // use Comm_split_type to split into shared memory subcommunicators
        MPI_Comm_split_type(this->mpi_comm, MPI_COMM_TYPE_SHARED, this->rank(), MPI_INFO_NULL, &o.mpi_comm);
#else
        // use `hostname` (dirty workaround) to split communicator
        char p_name[MPI_MAX_PROCESSOR_NAME+1];
        int p_len;
        MPI_Get_processor_name(p_name, &p_len);
        p_name[p_len] = '\0'; // make string NULL-terminated (if not yet so)
        std::string str(p_name);

        // hash hostname
        std::hash<std::string> hash_func;
        size_t hash = hash_func(str);
        int h = static_cast<int>(hash % std::numeric_limits<int>::max());

        // split communicator by hash
        comm hash_comm;
        MPI_Comm_split(this->mpi_comm, h, this->rank(), &hash_comm.mpi_comm);
	hash_comm.init_ranksize();

        // potentially further split in case of hash collisions
        // check if all names are the same per subcomm
        std::string next_name = mxx::left_shift(str, hash_comm);
        bool same = hash_comm.rank() == hash_comm.size()-1 || str == next_name;
        if (!mxx::all_of(same, hash_comm)) {
            // there's a hash collision => manually resolve
            // gather all strings to rank 0, process, and scatter assigned colors
            std::vector<size_t> recv_sizes = mxx::gather(static_cast<size_t>(p_len+1), 0, hash_comm);
            mxx::local_exscan_inplace(recv_sizes.begin(), recv_sizes.end());
            std::vector<char> hostnamesarr = mxx::gatherv(p_name, p_len+1, recv_sizes, 0, hash_comm);
            std::vector<int> colors;
            if (hash_comm.rank() == 0) {
                colors.resize(hash_comm.size());
                std::map<std::string, int> hostnames;
                int color = 0;
                for (int i = 0; i < hash_comm.size(); ++i) {
                    std::string n(&hostnamesarr[recv_sizes[i]]);
                    if (hostnames.find(n) == hostnames.end()) {
                        hostnames[n] = ++color;
                    }
                    colors[i] = hostnames[n];
                }
            }
            int my_color = mxx::scatter_one(colors, 0, hash_comm);
            // split by color
            MPI_Comm_split(hash_comm.mpi_comm, my_color, this->rank(), &o.mpi_comm);
        } else {
            o = std::move(hash_comm);
        }
#endif
        o.init_ranksize();
        o.do_free = true;
        return o;
}

} // namespace mxx

#endif // MXX_COMM_DEF_HPP

#endif
#endif
#endif
#endif
