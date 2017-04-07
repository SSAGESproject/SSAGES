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
 * @file    comm.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements a wrapper for MPI_Comm.
 */

#ifndef MXX_COMM_FWD_HPP
#define MXX_COMM_FWD_HPP

#include <mpi.h>

#include <functional>
#include <limits>

#include "common.hpp"
#include "future.hpp"
#include "datatypes.hpp"

namespace mxx {

static constexpr int any_tag = MPI_ANY_TAG;
static constexpr int any_source = MPI_ANY_SOURCE;


// forward declaration for streams
template <typename CharT, class Traits = std::char_traits<CharT> >
class sync_basic_ostream;

class comm {
public:
    /// Default constructor defaults to COMM_WORLD
    comm() : mpi_comm(MPI_COMM_WORLD) {
        MPI_Comm_size(mpi_comm, &m_size);
        MPI_Comm_rank(mpi_comm, &m_rank);
        do_free = false;
    }

    /// Taking an MPI_Comm object, user is responsible for freeing
    /// the communicator
    comm(const MPI_Comm& c) {
        mpi_comm = c;
        do_free = false;
        if (c != MPI_COMM_NULL)
            init_ranksize();
    }

    // disable copying
    comm(const comm& o) = delete;
    // disable copy assignment
    comm& operator=(const comm& o) = delete;

    /// Move constructor
    comm(comm&& o) : mpi_comm(o.mpi_comm), m_size(o.m_size), m_rank(o.m_rank), do_free(o.do_free) {
        o.mpi_comm = MPI_COMM_NULL;
        o.m_size = o.m_rank = 0;
        o.do_free = false;
    }

    /// Move assignment
    comm& operator=(comm&& o) {
        if (&o != this) {
            free();
            mpi_comm = o.mpi_comm;
            m_size = o.m_size;
            m_rank = o.m_rank;
            do_free = o.do_free;
            o.mpi_comm = MPI_COMM_NULL;
            o.m_size = o.m_rank = 0;
            o.do_free = false;
        }
        return *this;
    }

    /**
     * @brief   Explicity duplicate this communicator into a new communicator
     *          object. This is a collective call.
     *
     * @note    This is a collective operation and has to be called by all
     *          processes in this communicator.
     *
     * @return  The new, duplicated communicator.
     */
    comm copy() const {
        comm o;
        MPI_Comm_dup(mpi_comm, &o.mpi_comm);
        o.init_ranksize();
        o.do_free = true;
        return o;
    }

    /**
     * @brief   Splits the communicator into multiple sub-communicators, one for each color.
     *
     * @note    This is a collective operation and has to be called by all
     *          processes in this communicator.
     *
     * @return  The subcommunicator object.
     */
    comm split(int color) const {
        comm o;
        MPI_Comm_split(this->mpi_comm, color, this->rank(), &o.mpi_comm);
        o.init_ranksize();
        o.do_free = true;
        return o;
    }

    /**
     * @brief   Splits the communicator into multiple sub-communicators, one for each color,
     *          the order or ranks is assigned according to the given `key`.
     *
     * @note    This is a collective operation and has to be called by all
     *          processes in this communicator.
     *
     * @return  The subcommunicator object.
     */
    comm split(int color, int key) const {
        comm o;
        MPI_Comm_split(this->mpi_comm, color, key, &o.mpi_comm);
        o.init_ranksize();
        o.do_free = true;
        return o;
    }

    /**
     * @brief   Executes a given function with only a subset of processes.
     *
     * The given function should be of the signature: void(const mxx::comm&)
     * This communicator is split with the boolean condition.
     * The given function is called only for those processes for which the
     * boolean condition was true
     *
     * @tparam Func An object/function with operator()(const mxx::comm&).
     * @param cond  Ranks with `true` are executing the given function in
     *              a subcommunicator.
     * @param f     An object/function with operator()(const mxx::comm&), which
     *              gets called on ranks with `cond == true`.
     */
    template <typename Func>
    void with_subset(bool cond, Func f) const;

    template <typename Func>
    void exlusively_in_order_do(Func f) const {
        // executes f() in the order of ranks
        if (this->rank() > 0) {
            this->recv<int>(this->rank()-1, 345);
        }
        // call function
        f();
        // tell next processor to go ahead
        if (this->rank() != this->size()-1) {
            this->send(0, this->rank()+1, 345);
        }
    }

    // returns synchronized stream object for this communicator.
    // The stream object's destructor contains a collective operation
    // for synchronized cout/cerr output
    sync_basic_ostream<char> sync_cout() const;
    sync_basic_ostream<char> sync_cerr() const;


    /**
     * @brief   Returns a new communicator which is the reverse of this.
     *
     * @note    This is a collective operation and has to be called by all
     *          processes in this communicator.
     *
     * @return  The reverse communicator.
     */
    comm reverse() const {
        comm o;
        MPI_Comm_split(this->mpi_comm, 0, this->size() - this->rank(), &o.mpi_comm);
        o.init_ranksize();
        o.do_free = true;
        return o;
    }

    /**
     * @brief   Splits this communicator into subcommunicators, one for each
     *          node/shared memory accessible regions.
     *
     * @note    This is a collective operation and has to be called by all
     *          processes in this communicator.
     *
     * @return  The subcommunicator object.
     */
    comm split_shared() const;

public:

    /// Implicit conversion to MPI_Comm
    operator MPI_Comm() const {
        return mpi_comm;
    }

    /// Destructor (frees the MPI_Comm object)
    virtual ~comm() {
        free();
    }

    /// Returns the size of the communicator
    int size() const {
        return m_size;
    }

    /// Returns the rank of this process in the communicator
    int rank() const {
        return m_rank;
    }

    /// Collective barrier call for all processes in `this` communicator.
    void barrier() const {
        MPI_Barrier(this->mpi_comm);
    }

private:
    void free() {
        if (!is_builtin() && do_free) {
            MPI_Comm_free(&mpi_comm);
        }
    }

    bool is_builtin(const MPI_Comm& c) const {
        return c == MPI_COMM_WORLD
            || c == MPI_COMM_SELF
            || c == MPI_COMM_NULL;
    }

    bool is_builtin() const {
        return is_builtin(mpi_comm);
    }

    // initiate the rank and size based on the communicator
    void init_ranksize() {
        MPI_Comm_size(mpi_comm, &m_size);
        MPI_Comm_rank(mpi_comm, &m_rank);
    }

private:
    /// The MPI Communicator being wrapped
    MPI_Comm mpi_comm;
    /// The size of the communicator
    int m_size;
    /// This process' rank in this communicator
    int m_rank;
    /// Whether or not to free the MPI_Comm object upon destruction of `this`.
    bool do_free;

public:

    /**************
     *  MPI send  *
     **************/

    /// Send a basic (fixed size) datatype
    template <typename T>
    inline void send(const T& msg, int dest, int tag = 0) const {
        MXX_ASSERT(sizeof(T) < mxx::max_int);
        MXX_ASSERT(0 <= dest && dest < this->size());
        mxx::datatype dt = mxx::get_datatype<T>();
        MPI_Send(const_cast<T*>(&msg), 1, dt.type(), dest, tag, this->mpi_comm);
    }

    /// send a block of memory of a specified base datatype
    template <typename T>
    inline void send(const T* msg, size_t size, int dest, int tag = 0) const {
        MXX_ASSERT(0 <= dest && dest < this->size());
        MXX_ASSERT(msg != nullptr);
        if (size < mxx::max_int) {
            mxx::datatype dt = mxx::get_datatype<T>();
            MPI_Send(const_cast<T*>(msg), size, dt.type(), dest, tag, this->mpi_comm);
        } else {
            mxx::datatype dt = mxx::get_datatype<T>().contiguous(size);
            MPI_Send(const_cast<T*>(msg), 1, dt.type(), dest, tag, this->mpi_comm);
        }
    }

    /// Send a std::vector
    template <typename T>
    inline void send(const std::vector<T>& msg, int dest, int tag = 0) const {
        this->send(&msg[0], msg.size(), dest, tag);
    }

    /// Send a std::string
    template <typename CharT, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT> >
    inline void send(const std::basic_string<CharT, Traits, Alloc>& msg, int dest, int tag = 0) const {
        this->send(&msg[0], msg.size(), dest, tag);
    }

    /// pass string literal to comm::send
    template <typename CharT, size_t N>
    inline void send(const CharT(&msg)[N], int dest, int tag = 0) const {
        MXX_ASSERT(0 <= dest && dest < this->size());
        if (N < mxx::max_int) {
            mxx::datatype dt = mxx::get_datatype<CharT>();
            MPI_Send(const_cast<CharT*>(msg), N-1, dt.type(), dest, tag, this->mpi_comm);
        } else {
            mxx::datatype dt = mxx::get_datatype<CharT>().contiguous(N-1);
            MPI_Send(const_cast<CharT*>(msg), 1, dt.type(), dest, tag, this->mpi_comm);
        }
    }


    /***********************
     *  Non-blocking Send  *
     ***********************/

    /// Send a basic (fixed size) datatype
    template <typename T>
    inline mxx::future<void> isend(const T& msg, int dest, int tag = 0) const {
        MXX_ASSERT(sizeof(T) < mxx::max_int);
        MXX_ASSERT(0 <= dest && dest < this->size());
        mxx::datatype dt = mxx::get_datatype<T>;
        mxx::future_builder<void> f;
        MPI_Isend(const_cast<T*>(&msg), 1, dt.type(), dest, tag, this->mpi_comm, &f.add_request());
        return f.get_future();
    }

    /// send a block of memory of a specified base datatype
    template <typename T>
    inline mxx::future<void> isend(const T* msg, size_t size, int dest, int tag = 0) const {
        MXX_ASSERT(0 <= dest && dest < this->size());
        mxx::future_builder<void> f;
        if (size < mxx::max_int) {
            mxx::datatype dt = mxx::get_datatype<T>();
            MPI_Isend(const_cast<T*>(msg), size, dt.type(), dest, tag, this->mpi_comm, &f.add_request());
        } else {
            mxx::datatype dt = mxx::get_datatype<T>().contiguous(size);
            MPI_Isend(const_cast<T*>(msg), 1, dt.type(), dest, tag, this->mpi_comm, &f.add_request());
        }
        return f.get_future();
    }

    /// send a std::vector
    template <typename T>
    inline mxx::future<void> isend(const std::vector<T>& msg, int dest, int tag = 0) const {
        return this->isend(&msg[0], msg.size(), dest, tag);
    }

    /// Send a std::string
    template <typename CharT, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT> >
    inline mxx::future<void> isend(const std::basic_string<CharT,Traits,Alloc>& msg, int dest, int tag = 0) const {
        return this->isend(&msg[0], msg.size(), dest, tag);
    }


    /**********
     *  Recv  *
     **********/
    // TODO: recv should also return a `status`, so that we can known which tag and src the message is from

    // recv single element and return
    template <typename T>
    T recv(int src, int tag = mxx::any_tag) const;

    // recv into given buffer of size 1
    template <typename T>
    void recv_into(T& buffer, int src, int tag = mxx::any_tag) const;

    // recv into given buffer of given size
    template <typename T>
    void recv_into(T* buffer, size_t count, int src, int tag = mxx::any_tag) const;

    // recv into vector and return
    template <typename T, class Alloc = std::allocator<T>>
    std::vector<T, Alloc> recv_vec(size_t size, int src, int tag = mxx::any_tag) const;

    // recv into string and return
    template <typename CharT = char, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT> >
    std::basic_string<CharT, Traits, Alloc> recv_str(size_t size, int src, int tag = mxx::any_tag) const;


    /***************************
     *  Non-blocking receives  *
     ***************************/

    // recv single element and return
    template <typename T>
    mxx::future<T> irecv(int src, int tag = mxx::any_tag) const;

    // recv into given buffer of size 1
    // TODO: should this function be disabled for strict buffer protection?
    template <typename T>
    mxx::future<void> irecv_into(T& buffer, int src, int tag = mxx::any_tag) const;

    // recv into given buffer of given size
    template <typename T>
    mxx::future<void> irecv_into(T* buffer, size_t count, int src, int tag = mxx::any_tag) const;

    // recv into vector and return
    template <typename T, class Alloc = std::allocator<T>>
    mxx::future<std::vector<T, Alloc> > irecv_vec(size_t size, int src, int tag = mxx::any_tag) const;

    // recv into string and return
    template <typename CharT = char, class Traits = std::char_traits<CharT>, class Alloc = std::allocator<CharT> >
    mxx::future<std::basic_string<CharT, Traits, Alloc> > irecv_str(size_t size, int src, int tag = mxx::any_tag) const;
};

template <typename T>
struct recv_impl {
    static inline void do_recv_into(int src, int tag, T& buf, MPI_Comm c) {
        mxx::datatype dt = mxx::get_datatype<T>();
        MPI_Recv(&buf, 1, dt.type(), src, tag, c, MPI_STATUS_IGNORE);
    }
    static inline T do_recv(int src, int tag, MPI_Comm c) {
        T result;
        do_recv_into(src, tag, result, c);
        return result;
    }
};

// receive previously unknown sized contiguous container class (std::vector or std::string)
// uses:
//  - `Container.resize(size)`  to resize the container to the size received
//  - `&Container.operator[0]`  for getting the address of the first element
//  - `Container::value_type`   to get the underlying data type
//  The Container thus must be layed out contiguously in memory
template <typename Container>
struct recv_container_impl {
    static inline void do_recv_into(int src, int tag, Container& buf, MPI_Comm c) {
        // TODO: how do I do this in async?
        typedef typename Container::value_type T;
        mxx::datatype dt = mxx::get_datatype<T>();
        MPI_Status stat;
#if MPI_VERSION >= 3
        // threadsafe version with MProbe and MRecv (only if MPI-3)
        // and safe for > INT_MAX receives
        MPI_Message msg;
        MPI_Mprobe(src, tag, c, &msg, &stat);
        size_t size;
        MPI_Count count;
        MPI_Get_elements_x(&stat, dt.type(), &count);
        // TODO: (maybe) use a different mechanism for determining the
        //               number of basic elements per item of type `T`,
        //               e.g. by sending single element to MPI_COMM_SELF and MPI_Get_elements_x on status
        //               and do so the first time a type gets created
        //               and then cache the information in the datatype
        MXX_ASSERT(count % mxx::datatype_builder<T>::num_basic_elements() == 0);
        size = count / mxx::datatype_builder<T>::num_basic_elements();
        if (buf.size() != size)
            buf.resize(size);
        // receive into buffer
        if (size < mxx::max_int) {
            mxx::datatype dt = mxx::get_datatype<T>();
            MPI_Mrecv(const_cast<T*>(&buf[0]), size, dt.type(), &msg, &stat);
        } else {
            mxx::datatype dt = mxx::get_datatype<T>().contiguous(size);
            MPI_Mrecv(const_cast<T*>(&buf[0]), 1, dt.type(), &msg, &stat);
        }
#else
        // TODO: use communicator lock for thread safety??
        // TODO: not safe for messages with sizes larger than INT_MAX
        MPI_Probe(src, tag, c, &stat);
        int count;
        MPI_Get_count(&stat, dt.type(), &count);
        MXX_ASSERT(count >= 0); /* this assertion should fail if the message size was too large */
        if (buf.size() != count)
            buf.resize(count);
        MPI_Recv(&buf[0], count, dt.type(), stat.MPI_SOURCE, stat.MPI_TAG, c, MPI_STATUS_IGNORE);
#endif
    }
};

// template specialize for std::vector (variable sized!)
template <typename T, class Alloc>
struct recv_impl<std::vector<T, Alloc> > {
    static inline void do_recv_into(int src, int tag, std::vector<T>& buf, MPI_Comm c) {
        recv_container_impl<std::vector<T, Alloc> >::do_recv_into(src, tag, buf, c);
    }
    static inline std::vector<T, Alloc> do_recv(int src, int tag, MPI_Comm c) {
        std::vector<T, Alloc> result;
        recv_container_impl<std::vector<T, Alloc> >::do_recv_into(src, tag, result, c);
        return result;
    }
};
// template specialize for std::string (variable sized!)
template <typename CharT, class Traits, class Alloc>
struct recv_impl<std::basic_string<CharT, Traits, Alloc> > {
    static inline void do_recv_into(int src, int tag, std::basic_string<CharT, Traits, Alloc>& buf, MPI_Comm c) {
        recv_container_impl<std::basic_string<CharT, Traits, Alloc> >::do_recv_into(src, tag, buf, c);
    }
    static inline std::basic_string<CharT, Traits, Alloc> do_recv(int src, int tag, MPI_Comm c) {
        std::basic_string<CharT, Traits, Alloc> result;
        recv_container_impl<std::basic_string<CharT, Traits, Alloc> >::do_recv_into(src, tag, result, c);
        return result;
    }
};

// Regular templated receive returning the received data
// supports std::string and std::vector additionally to supported mxx::datatypes
// std::string and std::vector is only supported in blocking version
template <typename T>
inline T comm::recv(int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    return recv_impl<T>::do_recv(src, tag, this->mpi_comm);
}

template <typename T>
inline void comm::recv_into(T& buffer, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    return recv_impl<T>::do_recv_into(src, tag, buffer, this->mpi_comm);
}


template <typename T>
void comm::recv_into(T* buffer, size_t size, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    MXX_ASSERT(buffer != nullptr);
    if (size < mxx::max_int) {
        mxx::datatype dt = mxx::get_datatype<T>();
        MPI_Recv(const_cast<T*>(buffer), size, dt.type(), src, tag, this->mpi_comm, MPI_STATUS_IGNORE);
    } else {
        mxx::datatype dt = mxx::get_datatype<T>().contiguous(size);
        MPI_Recv(const_cast<T*>(buffer), 1, dt.type(), src, tag, this->mpi_comm, MPI_STATUS_IGNORE);
    }
}

// recv into vector and return
template <typename T, class Alloc>
inline std::vector<T, Alloc> comm::recv_vec(size_t size, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    std::vector<T, Alloc> result(size);
    this->recv_into(&result[0], size, src, tag);
    return result;
}

// recv into string and return
template <typename CharT, class Traits, class Alloc>
inline std::basic_string<CharT, Traits, Alloc> comm::recv_str(size_t size, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    std::basic_string<CharT, Traits, Alloc> result;
    result.resize(size);
    this->recv_into(&result[0], size, src, tag);
    return result;
}

/***********
 *  Irecv  *
 ***********/

template <typename T>
inline mxx::future<T> comm::irecv(int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    mxx::future_builder<T> f;
    mxx::datatype dt = mxx::get_datatype<T>();
    MPI_Irecv(f.data(), 1, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    return std::move(f.get_future());
}

// only for fixed size types, doesn't support std::vector or std::basic_string
// use `irecv_vec` or `irecv_str` instead
template <typename T>
inline mxx::future<void> comm::irecv_into(T& buffer, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    mxx::future_builder<void> f;
    mxx::datatype dt = mxx::get_datatype<T>();
    MPI_Irecv(&buffer, 1, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    return f.get_future();
}

// recv into given buffer of given size
template <typename T>
inline mxx::future<void> comm::irecv_into(T* buffer, size_t count, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    MXX_ASSERT(buffer != nullptr);
    mxx::future_builder<void> f;
    if (count < mxx::max_int) {
        mxx::datatype dt = mxx::get_datatype<T>();
        MPI_Irecv(buffer, count, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    } else {
        mxx::datatype dt = mxx::get_datatype<T>().contiguous(count);
        MPI_Irecv(buffer, 1, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    }
    return f.get_future();
}

// recv into vector and return
template <typename T, class Alloc>
inline mxx::future<std::vector<T, Alloc> > comm::irecv_vec(size_t size, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    mxx::future_builder<std::vector<T, Alloc> > f;
    f.data()->resize(size);
    if (size < mxx::max_int) {
        mxx::datatype dt = mxx::get_datatype<T>();
        MPI_Irecv(&(f.data()->front()), size, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    } else {
        mxx::datatype dt = mxx::get_datatype<T>().contiguous(size);
        MPI_Irecv(&(f.data()->front()), 1, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    }
    return f.get_future();
}

// recv into string and return
template <typename CharT, class Traits, class Alloc>
inline mxx::future<std::basic_string<CharT, Traits, Alloc> > comm::irecv_str(size_t size, int src, int tag) const {
    MXX_ASSERT(0 <= src && src < this->size());
    mxx::future_builder<std::basic_string<CharT, Traits, Alloc> > f;
    f.data()->resize(size);
    if (size < mxx::max_int) {
        mxx::datatype dt = mxx::get_datatype<CharT>();
        MPI_Irecv(&(f.data()->front()), size, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    } else {
        mxx::datatype dt = mxx::get_datatype<CharT>().contiguous(size);
        MPI_Irecv(&(f.data()->front()), 1, dt.type(), src, tag, this->mpi_comm, &f.add_request());
    }
    return std::move(f.get_future());
}

} // namespace mxx

#endif // MXX_COMM_HPP
