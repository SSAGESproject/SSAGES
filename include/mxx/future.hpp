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
 * @file    future.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements request wrappers and mxx::future.
 */

#ifndef MXX_FUTURE_HPP
#define MXX_FUTURE_HPP


// MPI include
#include <mpi.h>

#include <vector>
#include <memory>


namespace mxx
{

template <typename T>
class future;

template <typename T>
class future_builder;


class requests {
public:
    requests() : m_requests() {}
    requests(MPI_Request req) : m_requests(1, req) {}
    requests(const requests& req) = default;
    requests(requests&& req) = default;
    requests& operator=(const requests& req) = default;
    requests& operator=(requests&& req) = default;

    void append(MPI_Request req) {
        m_requests.push_back(req);
    }

    MPI_Request& add() {
        m_requests.push_back(MPI_Request());
        return m_requests.back();
    }

    MPI_Request& back() {
        return m_requests.back();
    }

    const MPI_Request& back() const {
        return m_requests.back();
    }

    MPI_Request& operator[](std::size_t i){
        return m_requests[i];
    }

    const MPI_Request& operator[](std::size_t i) const {
        return m_requests[i];
    }

    void waitall() {
        MPI_Waitall(m_requests.size(), &m_requests[0], MPI_STATUSES_IGNORE);
    }

    void wait() {
        this->waitall();
    }

    // adds all requests from `r` into this requests object
    void insert(requests& r) {
        // `steal` the requests
        m_requests.insert(m_requests.end(), r.m_requests.begin(), r.m_requests.end());
        r.m_requests.clear();
    }

    bool test() {
        int flag;
        MPI_Testall(m_requests.size(), &m_requests[0], &flag, MPI_STATUSES_IGNORE);
        return flag != 0;
    }

    // TODO: functions to access/return `MPI_Status`

    virtual ~requests() {
        this->waitall();
    }

private:
    std::vector<MPI_Request> m_requests;
};

/// Combines MPI request and received data storage similar to std::future
/// Calling .get() will first MPI_Wait and then std::move the data out of
/// the mxx::future
namespace impl {
template <typename T>
class future_base {
public:
    /// wrapped type
    typedef typename std::remove_reference<T>::type value_type;

    // disable copying
    future_base(const future_base& f) = delete;
    future_base& operator=(const future_base& f) = delete;

    // default move construction and assignment
    future_base(future_base&& f) = default;
    future_base& operator=(future_base&& f) = default;

    /// Returns `true` if the result is available
    bool valid() {
        if (!m_valid)
            m_valid = m_req.test();
        return m_valid;
    }

    /// blocks until the result becomes available
    void wait() {
        if (!m_valid)
            m_req.wait();
        m_valid = true;
        m_ever_valid = true;
    }

    virtual ~future_base() {
        // check if this has ever been valid. If not: wait(), otherwise
        // destruct the m_data member (happens anyway)
        if (!m_ever_valid && !m_valid)
            wait();
    }

protected:

    /// Default construction creates the output memory space (for MPI to write
    /// into).
    /// Only for friends!
    future_base() : m_valid(false), m_ever_valid(false) {}

    MPI_Request& add_request() {
        return m_req.add();
    }

    friend class mxx::future_builder<T>;

protected:
    bool m_valid;
    bool m_ever_valid;
    requests m_req;
};

}

template <typename T>
class future : public impl::future_base<T> {
public:
    /// wrapped type
    typedef typename std::remove_reference<T>::type value_type;

    // disable copying
    future(const future& f) = delete;
    future& operator=(const future& f) = delete;

    // default move construction and assignment
    future(future&& f) = default;
    future& operator=(future&& f) = default;

    value_type get() {
        this->wait();
        this->m_valid = false;
        return std::move(*m_data);
    }

protected:
    // functions only accessible by the future_builder
    value_type* data() {
        return m_data.get();
    }

    /// Default construction creates the output memory space (for MPI to write
    /// into).
    /// Only for friends!
    future() : impl::future_base<T>(), m_data(new T()) {}

    friend class mxx::future_builder<T>;

protected:
    typedef std::unique_ptr<value_type> ptr_type;
    ptr_type m_data;
};


// template specialization for <void>
template <>
class future<void> : public impl::future_base<void> {
public:
    /// wrapped type
    typedef void value_type;

    // disable copying
    future(const future& f) = delete;
    future& operator=(const future& f) = delete;

    // default move construction and assignment
    future(future&& f) = default;
    future& operator=(future&& f) = default;

    void get() {
        wait();
        this->m_valid = false;
    }

protected:
    // functions only accessible by the future_builder
    void* data() {
        return nullptr;
    }

    /// Default construction creates the output memory space (for MPI to write
    /// into).
    /// Only for friends!
    future() : impl::future_base<void>() {}

    // declare `builder` as friend
    friend class mxx::future_builder<void>;
};


// similar to a std::promise, this is a friend of mxx::future and
// enables mxx functions to build a mxx::future containing the
// yet to be written to buffers
template <typename T>
class future_builder {
public:
    typedef typename  mxx::future<T>::value_type value_type;

    future_builder() : m_valid(true), m_future() {}

    MPI_Request& add_request() {
        return m_future.add_request();
    }

    value_type* data() {
        return m_future.data();
    }

    mxx::future<T> get_future() {
        m_valid = false;
        return std::move(m_future);
    }

private:
    bool m_valid;
    mxx::future<T> m_future;
};


} // namespace mxx

#endif // MXX_FUTURE_HPP
