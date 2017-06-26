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
 * @file    benchmark.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements some MPI bandwidth benchmarks.
 */
#ifndef MXX_BENCHMARKS_HPP
#define MXX_BENCHMARKS_HPP

// MPI include
#include <mpi.h>


// C++ includes
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <cxx-prettyprint/prettyprint.hpp>

// MXX
#include "comm.hpp"
#include "datatypes.hpp"
#include "collective.hpp"
#include "timer.hpp"
#include "utils.hpp"

namespace mxx {

// a hierarchical communicator wrapper for MPI-MPI hybrid programming
// with 3 communicators: 1 global, 1 per node, 1 for node-masters
class hybrid_comm {
public:
    mxx::comm local;
    mxx::comm local_master;
    mxx::comm global;
    std::string node_name;

private:
    hybrid_comm()
        : local(MPI_COMM_NULL),
          local_master(MPI_COMM_NULL),
          global(MPI_COMM_NULL),
          node_name("") {}

public:
    hybrid_comm(const mxx::comm& c) :
        local(c.split_shared()),
        local_master(c.split(local.rank())),
        global(c.split(true, local_master.rank()*local.size() + local.rank())),
        node_name(mxx::get_processor_name()) {
        MXX_ASSERT(mxx::all_same(local.size()));
    }

    hybrid_comm split_by_node(int color) const {
        // split the processes but assert that each node is only in
        // one process
        MXX_ASSERT(mxx::all_same(color, local));
        hybrid_comm result;
        result.local = local.copy();
        result.global = global.split(color);
        result.local_master = local_master.split(color);
        return result;
    }

    // move constructor moves all members
    hybrid_comm(hybrid_comm&& o) = default;

    // executes only with a subset of nodes
    template <typename Func>
    void with_nodes(bool participate, Func func) const {
        int part = participate;
        MXX_ASSERT(mxx::all_same(part, local));
        if (mxx::all_same(part, global)) {
            if (participate) {
                func(*this);
            }
        } else {
            hybrid_comm hc(split_by_node(participate));
            if (participate) {
                func(hc);
            }
        }
        global.barrier();
    }

    int num_nodes() const {
        return global.size() / local.size();
    }

    int node_rank() const {
        MXX_ASSERT(global.rank() / local.size() == local_master.rank());
        return global.rank() / local.size();
    }

    bool is_local_master() const {
        return local.rank() == 0;
    }

    template <typename Func>
    void with_local_master(Func func) const {
        if (local.rank() == 0)
            func();
        global.barrier();
    }
};

// TODO: take vector/iterators instead of `n`
std::pair<double,double> bw_simplex(const mxx::comm& c, int partner, size_t n) {
    std::vector<size_t> vec(n);
    std::vector<size_t> result(n);
    std::generate(vec.begin(), vec.end(), std::rand);

    c.barrier();
    auto start = std::chrono::steady_clock::now();

    if (c.rank() < partner) {
        c.send(vec, partner);
    } else {
        c.recv_into(&result[0], vec.size(), partner);
    }
    auto end1 = std::chrono::steady_clock::now();
    if (c.rank() < partner) {
        c.recv_into(&result[0], vec.size(), partner);
    } else {
        c.send(vec, partner);
    }
    auto end2 = std::chrono::steady_clock::now();
    c.barrier();
    double time1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start).count();
    double time2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1).count();

    // calculate simplex bandwidth
    double bw1 = 8*(double)n*sizeof(size_t)/time1/1000.0;
    double bw2 = 8*(double)n*sizeof(size_t)/time2/1000.0;
    return std::pair<double,double>(bw1, bw2);
}

template <typename T>
double time_duplex(const mxx::comm& c, int partner, const std::vector<T>& sendvec, std::vector<T>& recvvec) {
    MXX_ASSERT(sendvec.size() == recvvec.size());

    size_t n = sendvec.size();
    c.barrier();
    auto start = std::chrono::steady_clock::now();
    mxx::datatype dt = mxx::get_datatype<size_t>();
    // sendrecv for full duplex
    MPI_Sendrecv(const_cast<T*>(&sendvec[0]), n, dt.type(), partner, 0, &recvvec[0], n, dt.type(), partner, 0, c, MPI_STATUS_IGNORE);
    //c.barrier();
    auto end = std::chrono::steady_clock::now();
    double time_p2p = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // calculate duplex bandwidth
    //double bw = 2*8*(double)n*sizeof(size_t)/time_p2p/1000.0;
    return time_p2p;
}

template <typename T>
double bw_duplex_per_node(const mxx::comm& c, int partner, const mxx::comm& smc, const std::vector<T>& sendvec, std::vector<T>& recvvec) {
    size_t n = sendvec.size();
    double time = time_duplex(c, partner, sendvec, recvvec);
    double maxtime_duplex = mxx::allreduce(time, mxx::max<double>(), smc);
    double bw = 2*8*n*sizeof(T)*smc.size()/maxtime_duplex/1000.0;
    return bw; // returns bandwidth in Gb/s
}


void bm(const mxx::comm& c, int partner, const mxx::comm& smc) {

    size_t n = 10000000;
    n /= smc.size();
    size_t MB = (n*sizeof(size_t))/1024/1024;
    std::vector<size_t> vec(n);
    std::vector<size_t> result(n);
    std::generate(vec.begin(), vec.end(), std::rand);

    std::string node_name = get_processor_name();

    double timed = time_duplex(c, partner, vec, result);
    double bw_send, bw_recv;
    std::tie(bw_send, bw_recv) = bw_simplex(c, partner, n);

    // calculate BW
    double maxtime_duplex = mxx::allreduce(timed, mxx::max<double>(), smc);
    double sum_bwd = 2*8*n*sizeof(size_t)*smc.size()/maxtime_duplex/1000.0;
    double sum_bw_send = mxx::allreduce(bw_send, std::plus<double>(), smc);
    double sum_bw_recv = mxx::allreduce(bw_recv, std::plus<double>(), smc);
    c.with_subset(smc.rank() == 0, [&](const mxx::comm& subcomm) {
        mxx::sync_cout(subcomm) << "[" << node_name << "]: Node BW Duplex = "
        << sum_bwd << " Gb/s (Simplex " << sum_bw_send << " Gb/s send, " << sum_bw_recv
        << " Gb/s recv) [" << MB*smc.size() << " MiB]" << std::endl;
    });
}

// returns a row of pairwise bw benchmark results on each process where smc.rank() == 0
std::vector<double> pairwise_bw_matrix(const hybrid_comm& hc) {
    int num_nodes = hc.num_nodes();
    size_t n = 1000000;
    std::vector<size_t> vec(n);
    std::vector<size_t> result(n);
    std::generate(vec.begin(), vec.end(), std::rand);
    n /= hc.local.size();
    // nodes get partnered
    int node_idx = hc.node_rank();
    std::vector<double> bw_row;
    if (hc.is_local_master())
        bw_row.resize(num_nodes);
    for (int dist = 1; dist < num_nodes; dist <<= 1) {
        if (hc.global.rank() == 0) {
            std::cout << "Benchmarking p2p duplex for dist = " << dist << std::endl;
        }
        int partner_block;
        if ((node_idx / dist) % 2 == 0) {
            // to left block
            partner_block = (node_idx/dist + 1)*dist;
        } else {
            partner_block = (node_idx/dist - 1)*dist;
        }
        int inblock_idx = node_idx % dist;
        for (int i = 0; i < dist; ++i) {
            int partner_node;
            if (partner_block >= node_idx)
                partner_node = partner_block + (inblock_idx + i) % dist;
            else
                partner_node = partner_block + (inblock_idx + (dist - i)) % dist;
            int partner = partner_node*hc.local.size() + hc.local.rank();
            // benchmark duplex with partner
            if (partner_node < num_nodes) {
                double bw = bw_duplex_per_node(hc.global, partner, hc.local, vec, result);
                if (hc.is_local_master()) {
                    bw_row[partner_node] = bw;
                }
            } else {
                // if this node doesn't participate in the benchmarking, it
                // still needs to call the barrier that is otherwise called
                // inside the time_duplex function
                hc.global.barrier();
            }
        }
    }
    return bw_row;
}

void print_bw_matrix_stats(const hybrid_comm& hc, const std::vector<double>& bw_row) {
    // print matrix
    hc.with_local_master([&](){
        // print matrix:
        mxx::sync_cout(hc.local_master) << "[" << hc.node_name << "]: " << std::fixed << std::setw(4) << std::setprecision(1) << std::setfill(' ') << bw_row <<  std::endl;

        // calculate avg bw per node
        // calc min and max
        double sum = 0.0;
        double max = 0.0;
        double min = std::numeric_limits<double>::max();
        int num_nodes = hc.num_nodes();
        int node_idx = hc.node_rank();
        for (int i = 0; i < num_nodes; ++i) {
            if (i != node_idx) {
                sum += bw_row[i];
                if (bw_row[i] > max)
                    max = bw_row[i];
                if (bw_row[i] < min)
                    min = bw_row[i];
            }
        }
        double global_max_bw = mxx::allreduce(max, mxx::max<double>(), hc.local_master);
        // output avg bw
        mxx::sync_cout(hc.local_master) << "[" << hc.node_name << "]: Average BW: " << sum / (num_nodes-1) <<  " Gb/s, Max BW: " << max << " Gb/s, Min BW: " << min << " Gb/s" << std::endl;
        // calc overall average
        double allsum = mxx::allreduce(sum, hc.local_master);
        if (hc.local_master.rank() == 0) {
            std::cout << "Overall Average BW: " << allsum / ((num_nodes-1)*(num_nodes-1)) << " Gb/s, Max: " << global_max_bw << " Gb/s" << std::endl;
        }
        // count how many connections are below 50% of max
        int count_bad = 0;
        int count_terrible = 0; // below 20%
        std::vector<int> bad_nodes;
        std::vector<int> terrible_nodes;
        for (int i = 0; i < num_nodes; ++i) {
            if (i != node_idx) {
                if (bw_row[i] < global_max_bw*0.5) {
                    ++count_bad;
                    bad_nodes.push_back(i);
                }
                if (bw_row[i] < global_max_bw*0.2) {
                    ++count_terrible;
                    terrible_nodes.push_back(i);
                }
            }
        }

        hc.local_master.with_subset(count_bad > 0, [&](const mxx::comm& outcomm) {
             mxx::sync_cout(outcomm) << "[" << hc.node_name << "]: " << count_bad << " bad connections (<50% max): " << bad_nodes << std::endl;
        });
        hc.local_master.with_subset(count_terrible > 0, [&](const mxx::comm& outcomm) {
             mxx::sync_cout(outcomm) << "[" << hc.node_name << "]: " << count_terrible << " terrible connections (<20% max): " << terrible_nodes << std::endl;
        });
    });
}


// all-pairwise bandwidth between all nodes
bool vote_off(const hybrid_comm& hc, int num_vote_off, const std::vector<double>& bw_row) {
    int num_nodes = hc.num_nodes();
    int node_idx = hc.node_rank();
    std::vector<bool> voted_off(num_nodes, false);
    hc.with_local_master([&](){
        // vote off bottlenecks: ie. nodes with lots of close to min connection
        // vote off slowest nodes
        // TODO: each node has `k` votes, vote off nodes with most votes
        //int k = hc.local_master.size()/2;
        double max = *std::max_element(bw_row.begin(), bw_row.end());
        double global_max_bw = mxx::allreduce(max, mxx::max<double>(), hc.local_master);
        for (int i = 0; i < num_vote_off; ++i) {
            // determine minimum of those not yet voted off
            double min = std::numeric_limits<double>::max();
            int off = 0;
            if (!voted_off[node_idx]) {
                for (int j = 0; j < num_nodes; ++j) {
                    if (j != node_idx && !voted_off[j]) {
                        if (bw_row[j] < min) {
                            min = bw_row[j];
                            off = j;
                        }
                    }
                }
            }
            double allmin = mxx::allreduce(min, mxx::min<double>(), hc.local_master);
            // vote only if min is within 0.1*max of allmin
            std::vector<int> vote_off(num_nodes, 0);
            std::vector<double> min_off(num_nodes, std::numeric_limits<double>::max());
            if (!voted_off[node_idx] && min <= allmin+0.1*global_max_bw) {
                vote_off[off] = 1;
                min_off[off] = min;
            }
            std::vector<int> total_votes = mxx::reduce(vote_off, 0, hc.local_master);
            std::vector<double> total_min = mxx::reduce(min_off, 0, mxx::min<double>(), hc.local_master);
            int voted_off_idx;
            if (hc.local_master.rank() == 0) {
                typedef std::tuple<int, int, double> vote_t;
                std::vector<vote_t> votes(num_nodes);
                for (int i = 0; i < num_nodes; ++i) {
                    votes[i] = vote_t(i, total_votes[i], total_min[i]);
                }
                // sort decreasing by number votes
                std::sort(votes.begin(), votes.end(), [](const vote_t& x, const vote_t& y) { return std::get<1>(x) > std::get<1>(y) || (std::get<1>(x) == std::get<1>(y) && std::get<2>(x) < std::get<2>(y));});

                // remove items with 0 votes
                auto it = votes.begin();
                while (std::get<1>(*it) > 0)
                    ++it;
                votes.resize(std::distance(votes.begin(), it));

                // print order
                std::cout << "Votes: " << votes << std::endl;
                std::cout << "Voting off node " << std::get<0>(votes[0]) << " with min " << std::get<2>(votes[0]) << "GiB/s of global min " << allmin << " GiB/s" << std::endl;

                voted_off_idx = std::get<0>(votes[0]);
            }
            MPI_Bcast(&voted_off_idx, 1, MPI_INT, 0, hc.local_master);
            voted_off[voted_off_idx] = true;
            // TODO: run all2all benchmark after each voting
        }
    });
    int votedoff = voted_off[node_idx];
    MPI_Bcast(&votedoff, 1, MPI_INT, 0, hc.local);
    return votedoff == 0;
}

inline mxx::sync_ostream sync_ofstream(const mxx::comm& comm, std::ofstream& os, int root = 0) {
    return comm.rank() == root ? sync_ostream(comm, root, os) : sync_ostream(comm, root);
}

void write_new_nodefile(const hybrid_comm& hc, bool participate, const std::string& filename) {
    hc.with_nodes(participate, [&](const hybrid_comm& subhc) {
        subhc.with_local_master([&](){
            // on rank 0 open filename as stream and write out via sync stream?
            std::ofstream ofile;
            if (subhc.global.rank() == 0) {
                ofile.open(filename);
            }
            mxx::sync_ofstream(subhc.local_master, ofile) << hc.node_name << std::endl;
            ofile.close();
        });
    });
}


void bw_all2all(const mxx::comm& c, const mxx::comm& smc) {
    // message size per target processor
    for (int k = 1; k <= 16; k <<= 1) {
        int m = k*1024;
        std::vector<size_t> els(m*c.size());
        std::generate(els.begin(), els.end(), std::rand);
        std::vector<size_t> rcv(m*c.size());
        mxx::datatype dt = mxx::get_datatype<size_t>();
        c.barrier();
        auto start = std::chrono::steady_clock::now();
        MPI_Alltoall(&els[0], m, dt.type(), &rcv[0], m, dt.type(), c);
        auto end = std::chrono::steady_clock::now();
        // time in microseconds
        double time_all2all = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        double max_time = mxx::allreduce(time_all2all, mxx::max<double>(), c);
        double min_time = mxx::allreduce(time_all2all, mxx::min<double>(), c);
        size_t bits_sendrecv = 2*8*sizeof(size_t)*m*(c.size() - smc.size());
        // bandwidth in Gb/s
        double bw = bits_sendrecv / max_time / 1000.0;
        if (c.rank() == 0) {
            std::cout << "All2all bandwidth: " << bw << " Gb/s [min=" << min_time/1000.0 << " ms, max=" << max_time/1000.0 << " ms, local_size=" << bits_sendrecv/1024/1024 << " MiB]" << std::endl;
        }
    }
}

void bw_all2all_char(const mxx::comm& c, const mxx::comm& smc) {
    // message size per target processor
    for (int k = 8; k <= 64; k <<= 1) {
        int m = k*1024;
        std::srand(13*c.rank());
        std::vector<size_t> send_counts(c.size());
        for (int i = 0; i < c.size(); ++i) {
            send_counts[i] = m;
        }
        size_t n = std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0));
        std::vector<char> els(n);
        std::generate(els.begin(), els.end(), std::rand);
        std::vector<size_t> recv_counts = mxx::all2all(send_counts, c);
        size_t recv_n = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));
        std::vector<char> rcv(recv_n);
        mxx::datatype dt = mxx::get_datatype<char>();
        c.barrier();
        auto start = std::chrono::steady_clock::now();
        mxx::all2allv(&els[0], send_counts, &rcv[0], recv_counts, c);
        //MPI_Alltoall(&els[0], m, dt.type(), &rcv[0], m, dt.type(), c);
        auto end = std::chrono::steady_clock::now();
        // time in microseconds
        double time_all2all = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        double max_time = mxx::allreduce(time_all2all, mxx::max<double>(), c);
        double min_time = mxx::allreduce(time_all2all, mxx::min<double>(), c);
        size_t bits_sendrecv = 2*8*sizeof(char)*m*(c.size() - smc.size());
        // bandwidth in Gb/s
        double bw = bits_sendrecv / max_time / 1000.0;
        if (c.rank() == 0) {
            std::cout << "All2all char bandwidth: " << bw << " Gb/s [min=" << min_time/1000.0 << " ms, max=" << max_time/1000.0 << " ms, local_size=" << bits_sendrecv/1024/1024 << " MiB]" << std::endl;
        }
    }
}

void bw_all2all_unaligned_char(const mxx::comm& c, const mxx::comm& smc, bool realign) {
    // message size per target processor
    for (int k = 1; k <= 128; k <<= 1) {
        int m = k*1024;
        std::srand(13*c.rank());

        // send counts
        std::vector<size_t> send_counts(c.size());
        for (int i = 0; i < c.size(); ++i) {
            send_counts[i] = m + std::rand() % 8;
        }
        size_t n = std::accumulate(send_counts.begin(), send_counts.end(), static_cast<size_t>(0));
        // generate input
        std::vector<char> els(n);
        std::generate(els.begin(), els.end(), std::rand);
        // original recv counts
        std::vector<size_t> recv_counts = mxx::all2all(send_counts, c);
        size_t recv_n = std::accumulate(recv_counts.begin(), recv_counts.end(), static_cast<size_t>(0));
        std::vector<char> rcv(recv_n);

        // time all2all
        c.barrier();
        auto start = std::chrono::steady_clock::now();
        /*
        if (realign) {
            char_all2allv(&els[0], send_counts, &rcv[0], recv_counts, c);
        } else {
        */
            mxx::all2allv(&els[0], send_counts, &rcv[0], recv_counts, c);
        //}
        //mxx::datatype dt = mxx::get_datatype<char>();
        //MPI_Alltoall(&els[0], m, dt.type(), &rcv[0], m, dt.type(), c);
        auto end = std::chrono::steady_clock::now();
        if (realign) {
            std::vector<char> rcv2(recv_n);
            mxx::all2allv(&els[0], send_counts, &rcv2[0], recv_counts, c);
            if (!(rcv == rcv2)) {
                std::cout << "[ERROR] Vectors are not same" << std::endl;
                std::cout << "rcv=" << rcv << std::endl;
                std::cout << "rcv2=" << rcv2 << std::endl;
            }
        }
        // time in microseconds
        double time_all2all = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        double max_time = mxx::allreduce(time_all2all, mxx::max<double>(), c);
        double min_time = mxx::allreduce(time_all2all, mxx::min<double>(), c);
        size_t bits_sendrecv = 2*8*sizeof(char)*m*(c.size() - smc.size());
        // bandwidth in Gb/s
        double bw = bits_sendrecv / max_time / 1000.0;
        if (c.rank() == 0) {
            std::cout << "All2all UNaligned char bandwidth: " << bw << " Gb/s [min=" << min_time/1000.0 << " ms, max=" << max_time/1000.0 << " ms, local_size=" << bits_sendrecv/1024/1024 << " MiB]" << std::endl;
        }
    }
}

void benchmark_nodes_bw_p2p(const mxx::comm& comm = mxx::comm()) {
    // pair with another node
    hybrid_comm hc(comm);
    int proc_per_node = hc.local.size();

    // assert same number processors per node
    if (!mxx::all_same(proc_per_node, comm)) {
        std::cerr << "Error: this benchmark assumes the same number of processors per node" << std::endl;
        MPI_Abort(comm, -1);
    }

    int num_nodes = hc.num_nodes();

    if (num_nodes % 2 != 0) {
        std::cerr << "Error: this benchmark assumes an even number of nodes" << std::endl;
        MPI_Abort(comm, -1);
    }

    //for (int local_p = 1; local_p <= proc_per_node; local_p <<= 1) {
    /*
    int local_p = sm_comm.size();
        bool participate =  true; //sm_comm.rank() < local_p;
        mxx::comm c = comm.split(participate, node_idx*local_p + sm_comm.rank());
        */
     //   mxx::comm smc = sm_comm.split(participate);

        if (true) {
            std::vector<double> bw_row = pairwise_bw_matrix(hc);
            print_bw_matrix_stats(hc, bw_row);
            bool part = vote_off(hc, 4, bw_row); // TODO: process result
            if (hc.global.rank() == 0)
                std::cout << "Before vote off: " << std::endl;
            bw_all2all(hc.global, hc.local);
            if (hc.global.rank() == 0)
                std::cout << "After vote off: " << std::endl;
            hc.with_nodes(part, [&](const hybrid_comm& subhc) {
                bw_all2all(subhc.global, subhc.local);
                bw_all2all_char(subhc.global, subhc.local);
                bw_all2all_unaligned_char(subhc.global, subhc.local, false);
                if (subhc.global.rank() == 0)
                    std::cout << "== With re-alignment" << std::endl;
                bw_all2all_unaligned_char(subhc.global, subhc.local, true);
            });
            write_new_nodefile(hc, part, "blah.nodes");

            /*
            if (c.rank() == 0) {
                std::cout << "Running with " << local_p << "/" << proc_per_node << " processes per node" << std::endl;
            }
            MXX_ASSERT(c.size() % 2 == 0);
            if (local_p > 1) {
                // intranode BW test
                if (c.rank() == 0)
                    std::cout << "Intranode BW test" << std::endl;
                int partner = (c.rank() % 2 == 0) ? c.rank() + 1 : c.rank() - 1;
                bm(c, partner, smc);
            }
            // 1) closest neighbor
            if (c.rank() == 0)
                std::cout << "Closest Neighbor BW test" << std::endl;
            int partner = (node_idx % 2 == 0) ? c.rank() + local_p : c.rank() - local_p;
            bm(c, partner, smc);

            // 2) furthest neighbor
            if (c.rank() == 0)
                std::cout << "Furthest Neighbor BW test" << std::endl;
            partner = (c.rank() < c.size()/2) ? c.rank() + c.size()/2 : c.rank() - c.size()/2;
            bm(c, partner, smc);
            */
        }
    //}
    // wait for other processes to finish the benchmarking
    //comm.barrier();
}

} // namespace mxx

#endif // MXX_BENCHMARKS_HPP
