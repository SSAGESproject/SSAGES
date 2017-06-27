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
 * @file    shift.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   MPI shift communication patterns (exchange of boundary elements).
 *
 * TODO:
 * - [ ] shifting of strings
 * - [ ] support shifting in more than 1 dimension using cartesian communicators
 */

#ifndef MXX_SHIFT_HPP
#define MXX_SHIFT_HPP

// MPI include
#include <mpi.h>

// C++ includes
#include <vector>
#include <memory>

// mxx includes
#include "datatypes.hpp"
#include "comm_fwd.hpp"
#include "future.hpp"


namespace mxx
{

template <typename T>
T right_shift(const T& t, const mxx::comm& comm = mxx::comm()) {
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 13;

    T left_value = T();
    MPI_Request recv_req;
    // if not last processor
    if (comm.rank() > 0) {
        MPI_Irecv(&left_value, 1, dt.type(), comm.rank()-1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() < comm.size()-1) {
        // send my most right element to the right
        MPI_Send(const_cast<T*>(&t), 1, dt.type(), comm.rank()+1, tag, comm);
    }
    if (comm.rank() > 0) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
    return left_value;
}

template <typename T>
void right_shift(const T* in, size_t n, T* out, const mxx::comm& comm = mxx::comm())
{
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 13;

    MPI_Request recv_req;
    // if not last processor
    if (comm.rank() > 0) {
        MPI_Irecv(out, n, dt.type(), comm.rank()-1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() < comm.size()-1) {
        // send my most right element to the right
        MPI_Send(const_cast<T*>(in), n, dt.type(), comm.rank()+1, tag, comm);
    }
    if (comm.rank() > 0) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
}

template <typename T>
std::vector<T> right_shift(const std::vector<T>& v, const mxx::comm& comm = mxx::comm())
{
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 13;
    // receive the size first
    std::vector<T> result;
    size_t left_size = right_shift(v.size(), comm);

    MPI_Request recv_req;
    // if not last processor
    if (comm.rank() > 0 && left_size > 0) {
        result.resize(left_size);
        MPI_Irecv(&result[0], left_size, dt.type(), comm.rank()-1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() < comm.size()-1 && v.size() > 0) {
        // send my most right element to the right
        MPI_Send(const_cast<T*>(&v[0]),v.size() , dt.type(), comm.rank()+1, tag, comm);
    }
    if (comm.rank() > 0 && left_size > 0) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
    return result;
}


template <typename T>
T left_shift(const T& t, const mxx::comm& comm = mxx::comm()) {
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 15;

    T right_value = T();
    T left_val = t;
    MPI_Request recv_req;
    // if not last processor
    if (comm.rank() < comm.size()-1) {
        MPI_Irecv(&right_value, 1, dt.type(), comm.rank()+1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() > 0) {
        // send my most right element to the right
        MPI_Send(const_cast<T*>(&left_val), 1, dt.type(), comm.rank()-1, tag, comm);
    }
    if (comm.rank() < comm.size()-1) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
    return right_value;
}

template <typename T>
void left_shift(const T* in, size_t n, T* out, const mxx::comm& comm = mxx::comm())
{
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 15;

    MPI_Request recv_req;
    // if not last processor
    if (comm.rank() < comm.size()-1) {
        MPI_Irecv(out, n, dt.type(), comm.rank()+1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() > 0) {
        // send my most right element to the right
        MPI_Send(const_cast<T*>(in), n, dt.type(), comm.rank()-1, tag, comm);
    }
    if (comm.rank() < comm.size()-1) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
}

// specialization for std::vector
template <typename T>
std::vector<T> left_shift(const std::vector<T>& v, const mxx::comm& comm = mxx::comm())
{
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 15;
    // receive the size first
    std::vector<T> result;
    size_t right_size = left_shift(v.size(), comm);

    MPI_Request recv_req;
    // if not last processor
    // TODO: replace with comm.send/ comm.recv which automatically will resolve
    // to BIG MPI if message size is too large
    if (comm.rank() < comm.size()-1 && right_size > 0) {
        result.resize(right_size);
        MPI_Irecv(&result[0], right_size, dt.type(), comm.rank()+1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() > 0 && v.size() > 0) {
        // send my most right element to the right
        MPI_Send(const_cast<T*>(&v[0]), v.size(), dt.type(), comm.rank()-1, tag, comm);
    }
    if (comm.rank() < comm.size()-1 && right_size > 0) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
    return result;
}

// specialization for std::string
template <typename CharT>
std::basic_string<CharT> left_shift(const std::basic_string<CharT>& str, const mxx::comm& comm = mxx::comm())
{
    // get datatype
    datatype dt = get_datatype<CharT>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 15;
    // receive the size first
    std::basic_string<CharT> result;
    size_t str_len = str.size();
    size_t right_size = left_shift(str_len, comm);

    MPI_Request recv_req;
    // if not last processor
    // TODO: replace with comm.send/ comm.recv which automatically will resolve
    // to BIG MPI if message size is too large
    if (comm.rank() < comm.size()-1 && right_size > 0) {
        result.resize(right_size);
        MPI_Irecv(&result[0], right_size, dt.type(), comm.rank()+1, tag,
                  comm, &recv_req);
    }
    // if not first processor
    if (comm.rank() > 0 && str.size() > 0) {
        // send my most right element to the right
        MPI_Send(const_cast<CharT*>(&str[0]), str.size(), dt.type(), comm.rank()-1, tag, comm);
    }
    if (comm.rank() < comm.size()-1 && right_size > 0) {
        // wait for the async receive to finish
        MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
    }
    return result;
}


template <typename T>
mxx::future<T> async_right_shift(const T& x, const mxx::comm& comm = mxx::comm()) {
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 15;

    mxx::future_builder<T> f;
    // if not first processor
    if (comm.rank() > 0) {
        MPI_Irecv(f.data(), 1, dt.type(), comm.rank()-1, tag,
                  comm, &f.add_request());
    }
    // if not last processor
    if (comm.rank() < comm.size()-1){
        // send my most right element to the right
        MPI_Isend(const_cast<T*>(&x), 1, dt.type(), comm.rank()+1,
                  tag, comm, &f.add_request());
    }

    return std::move(f.get_future());
}

template <typename T>
mxx::future<T> async_left_shift(const T& x, const mxx::comm& comm = mxx::comm()) {
    // get datatype
    datatype dt = get_datatype<T>();

    // TODO: handle tags with MXX (get unique tag function)
    int tag = 15;

    mxx::future_builder<T> f;
    // if not last processor
    if (comm.rank() < comm.size()-1) {
        MPI_Irecv(f.data(), 1, dt.type(), comm.rank()+1, tag,
                  comm, &f.add_request());
    }
    // if not first processor
    if (comm.rank() > 0) {
        // send my most right element to the right
        MPI_Isend(const_cast<T*>(&x), 1, dt.type(), comm.rank()-1,
                  tag, comm, &f.add_request());
    }
    return std::move(f.get_future());
}

} // namespace mxx


#define MXX_SHIFT_DONE
#include "comm_def.hpp"

#endif // MXX_SHIFT_HPP
