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
 * @file    reduction.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Reduction operations.
 */

#ifndef MXX_REDUCTION_HPP
#define MXX_REDUCTION_HPP

#include <mpi.h>

#include <vector>
#include <iterator>
#include <limits>
#include <functional>
#include <mutex>
#include <atomic>
#include <assert.h> // TODO: replace with own assert (calling MPI_Abort)
#include <cstring> // for memcpy

// mxx includes
#include "datatypes.hpp"
#include "comm_fwd.hpp"
#include "shift.hpp"

namespace mxx {

/*
 * Add any (key,value) pair to any MPI_Datatype using MPI caching (keyval)
 * functionality (used internally to mxx for adding std::function objects to
 * MPI datatypes, required for custom operators)
 */
template <typename K, typename T>
class attr_map {
private:
    // static "variables" wrapped into static member functions
    static int& key(){ static int k=0; return k; }
    static std::mutex& mut(){ static std::mutex m; return m; }
    static std::map<K, int>& keymap(){ static std::map<K, int> m; return m; }

    static int copy_attr(MPI_Datatype, int, void*, void *attribute_val_in, void *attribute_val_out, int *flag) {
        T** out = (T**) attribute_val_out;
        // copy via copy constructor
        *out = new T(*(T*)(attribute_val_in));
        *flag = 1;
        return MPI_SUCCESS;
    }

    static int del_attr(MPI_Datatype, int, void *attribute_val, void *) {
        delete (T*)attribute_val;
        return MPI_SUCCESS;
    }
public:
    static void set(const MPI_Datatype& dt, const K& k, const T& value) {
        std::lock_guard<std::mutex> lock(mut());
        if (keymap().find(k) == keymap().end()) {
            // create new keyval with MPI
            keymap()[k] = 0;
            MPI_Type_create_keyval(&attr_map<K,T>::copy_attr,
                                   &attr_map<K,T>::del_attr,
                                   &(keymap()[k]), (void*)NULL);
        }
        // insert new key into keymap
        // set the keyval pair
        T* val = new T(value); // copy construct
        MPI_Type_set_attr(dt, keymap()[k], (void*)val);
    }
    static T& get(const MPI_Datatype& dt, const K& k) {
        std::lock_guard<std::mutex> lock(mut());
        if (keymap().find(k) == keymap().end()) {
            throw std::out_of_range("Key not mapped.");
        }
        int key = keymap()[k];
        T* result;
        int flag;
        MPI_Type_get_attr(dt, key, &result, &flag);
        if (!flag) {
            throw std::out_of_range("Key not mapped.");
        }
        return *result;
    }
};

/*********************************************************************
 *                     User supplied functions:                      *
 *********************************************************************/

//  std::min and std::max are overloaded functions, and thus can't easily
//  be passed as functors.
//  Thus we define functors mxx::min<T> and mxx::max<T> with similar
//  declarations as std::plus, std::multiply etc

template <typename T>
struct min {
    inline const T& operator()(const T& x, const T& y) const {
        return std::min(x, y);
    }
};

template <typename T>
struct max {
    inline const T& operator()(const T& x, const T& y) const {
        return std::max(x, y);
    }
};

// an Op is not built-in in general:
template <typename T, typename Func>
struct get_builtin_op {
    static MPI_Op op(Func) {
        return MPI_OP_NULL;
    }
};

// define all C++ functors that map to a MPI builtin MPI_Op to directly map
// to that builtin MPI_Op
#define MXX_BUILTIN_OP(mpi_op, cpp_functor)                                    \
template <typename T>                                                          \
struct get_builtin_op<T, cpp_functor<T> > {                                    \
    static MPI_Op op(cpp_functor<T>) {                                         \
        return mpi_op;                                                         \
    }                                                                          \
};                                                                             \

MXX_BUILTIN_OP(MPI_SUM, std::plus);
MXX_BUILTIN_OP(MPI_PROD, std::multiplies);
MXX_BUILTIN_OP(MPI_LAND, std::logical_and);
MXX_BUILTIN_OP(MPI_LOR, std::logical_or);
MXX_BUILTIN_OP(MPI_BOR, std::bit_or);
MXX_BUILTIN_OP(MPI_BXOR, std::bit_xor);
MXX_BUILTIN_OP(MPI_BAND, std::bit_and);
MXX_BUILTIN_OP(MPI_MAX, mxx::max);
MXX_BUILTIN_OP(MPI_MIN, mxx::min);

#undef MXX_BUILTIN_OP


// for std::min/std::max functions
template <typename T>
struct get_builtin_op<T, const T&(*) (const T&, const T&)> {
    static MPI_Op op(const T& (*t)(const T&, const T&)) {
        // check if function is std::min or std::max
        if (t == static_cast<const T&(*)(const T&, const T&)>(std::min<T>)){
            return MPI_MIN;
        } else if (t == static_cast<const T&(*)(const T&, const T&)>(std::max<T>)){
            return MPI_MAX;
        } else {
            // otherwise return NULL
            return MPI_OP_NULL;
        }
    }
};

/**
 * @brief   Wrapps a binary combination/reduction operator for MPI use in
 *          custom operators.
 *
 * @note    This assumes that the operator is commutative.
 *
 * @tparam T    The input and ouput datatype of the binary operator.
 * @tparam IsCommutative    Whether or not the operation is commutative (default = true).
 */
template <typename T, bool IsCommutative = true>
class custom_op {
public:

    /**
     * @brief Creates a custom operator given a functor and the associated
     *        `MPI_Datatype`.
     *
     * @tparam Func     Type of the functor, can be a function pointer, lambda
     *                  function, or std::function or any object with a
     *                  `T operator(T& x, T& y)` member.
     * @param func      The instance of the functor.
     */
    template <typename Func>
    custom_op(Func func) : m_builtin(false) {
        if (mxx::is_builtin_type<T>::value) {
            // check if the operator is MPI built-in (in case the type
            // is also a MPI built-in type)
            MPI_Op op = get_builtin_op<T, Func>::op(std::forward<Func>(func));
            if (op != MPI_OP_NULL) {
                // this op is builtin, save it as such and don't copy built-in type
                m_builtin = true;
                m_op = op;
                mxx::datatype dt = mxx::get_datatype<T>();
                m_type_copy = dt.type();
            }
        }
        if (!m_builtin) {
            // create user function
            using namespace std::placeholders;
            m_user_func = std::bind(custom_op::custom_function<Func>,
                                  std::forward<Func>(func), _1, _2, _3, _4);
            // get datatype associated with the type `T`
            mxx::datatype dt = mxx::get_datatype<T>();
            // attach function to a copy of the datatype
            MPI_Type_dup(dt.type(), &m_type_copy);
            attr_map<int, func_t>::set(m_type_copy, 1347, m_user_func);
            // create op
            MPI_Op_create(&custom_op::mpi_user_function, IsCommutative, &m_op);
        }
    }


    /**
     * @brief   Returns the MPI_Datatype which has to be used in conjuction
     *          with the MPI_Op operator.
     *
     * The custom operators are wrapped into `std::function` objects and
     * saved/attached to a duplicated MPI_Datatype as MPI attribute.
     *
     * When MPI calls the custom user function, the MPI_Datatype is supplied
     * and thus the `std::function` object can be accessed and executed.
     *
     * @returns The `MPI_Datatype` which has to be used in conjuction with the
     *          `MPI_Op` returned by `get_op()` for all MPI reduction operations.
     */
    MPI_Datatype get_type() const {
        return m_type_copy;
    }

    /**
     * @brief   Returns the `MPI_Op` operator for reduction operations.
     *
     * @note
     * The MPI operator `MPI_Op` returned by this function can only be used in
     * conjunction with the `MPI_Datatype` returned by `get_type()`.
     *
     * @returns     The MPI operator as `MPI_Op` object.
     */
    MPI_Op get_op() const {
        return m_op;
    }

    /// Destructor: cleanup MPI objects
    virtual ~custom_op() {
        if (!m_builtin) {
            // clean-up (only if this wasn't a built-in MPI_Op)
            MPI_Op_free(&m_op);
            MPI_Type_free(&m_type_copy);
        }
    }
private:
    // Apply the user provided function to all elements passed by MPI.
    // The user provided function (lambda, function pointer, functor)
    // is bound to this function via std::bind, and the resulting object
    // saved in the MPI_Datatype
    template <typename Func>
    static void custom_function(Func func, void* invec, void* inoutvec, int* len, MPI_Datatype* dt) {
        if (*len > 1) {
            T* in = (T*) invec;
            T* inout = (T*) inoutvec;
            for (int i = 0; i < *len-1; ++i) {
                inout[i] = func(in[i], inout[i]);
            }
        }
        // only read and write the `true_extent` of the datatype
        // for the last item (otherwise we might encounter a memory error)
        // [see github OpenMPI issue #1462]
        T in_buf;
        T inout_buf;
        MPI_Aint true_extent, true_lb;
        MPI_Type_get_true_extent(*dt, &true_lb, &true_extent);
        std::memcpy((char*)&in_buf+true_lb, (char*)invec+(sizeof(T)*(*len-1)), true_extent);
        std::memcpy((char*)&inout_buf+true_lb, (char*)inoutvec+(sizeof(T)*(*len-1)), true_extent);
        // now do the operation on our local buffers
        inout_buf = func(in_buf, inout_buf);
        // copy the results back, again only to true_extent
        std::memcpy((char*)inoutvec+(sizeof(T)*(*len - 1)), (char*)&inout_buf+true_lb, true_extent);
    }
    // MPI custom Op function: (of type MPI_User_function)
    // This function is called from within MPI
    static void mpi_user_function(void* in, void* inout, int* n, MPI_Datatype* dt) {
        // get the std::function from the MPI_Datatype and call it
        func_t f = attr_map<int, func_t>::get(*dt, 1347);
        f(in, inout, n, dt);
    }

    // the std::function user function wrapper, which is called from the mpi user function
    typedef std::function<void(void*,void*,int*, MPI_Datatype*)> func_t;
    func_t m_user_func;
    /// Whether the MPI_Op is a builtin operator (e.g. MPI_SUM)
    bool m_builtin;
    /// The copy (Type_dup) of the MPI_Datatype to work on
    MPI_Datatype m_type_copy;
    /// The MPI user operator
    MPI_Op m_op;
};



/*********************************************************************
 *                Reductions                                         *
 *********************************************************************/


/*********************************************************************
 *                              Reduce                               *
 *********************************************************************/

template <typename T, typename Func>
inline void reduce(const T* in, size_t n, T* out, int root, Func func, const mxx::comm& comm = mxx::comm()) {
    // get custom op
    mxx::custom_op<T> op(std::forward<Func>(func));
    MPI_Reduce(const_cast<T*>(in), out, n, op.get_type(), op.get_op(), root, comm);
}

template <typename T>
inline void reduce(const T* in, size_t n, T* out, int root, const mxx::comm& comm = mxx::comm()) {
    reduce(in, n, out, root, std::plus<T>(), comm);
}

template <typename T, typename Func>
inline std::vector<T> reduce(const T* in, size_t n, int root, Func func, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result;
    if (comm.rank() == root)
        result.resize(n);
    reduce(in, n, &result[0], root, func, comm);
    return result;
}

template <typename T>
inline std::vector<T> reduce(const T* in, size_t n, int root, const mxx::comm& comm = mxx::comm()) {
    return reduce(in, n, root, std::plus<T>(), comm);
}

template <typename T, typename Func>
inline std::vector<T> reduce(const std::vector<T>& in, int root, Func func, const mxx::comm& comm = mxx::comm()) {
    return reduce(&in[0], in.size(), root, func, comm);
}

template <typename T>
inline std::vector<T> reduce(const std::vector<T>& in, int root, const mxx::comm& comm = mxx::comm()) {
    return reduce(in, root, std::plus<T>(), comm);
}

template <typename T, typename Func>
inline T reduce(const T& x, int root, Func func, const mxx::comm& comm = mxx::comm()) {
    // get custom op (and type for custom op)
    mxx::custom_op<T> op(std::forward<Func>(func));
    T result = T();
    MPI_Reduce(const_cast<T*>(&x), &result, 1, op.get_type(), op.get_op(), root, comm);
    return result;
}

template <typename T>
inline T reduce(const T& x, int root, const mxx::comm& comm = mxx::comm()) {
    return reduce(x, root, std::plus<T>(), comm);
}


/*********************************************************************
 *                             Allreduce                             *
 *********************************************************************/

template <typename T, typename Func>
inline void allreduce(const T* in, size_t n, T* out, Func func, const mxx::comm& comm = mxx::comm()) {
    // get custom op
    mxx::custom_op<T> op(std::forward<Func>(func));
    MPI_Allreduce(const_cast<T*>(in), out, n, op.get_type(), op.get_op(), comm);
}

template <typename T>
inline void allreduce(const T* in, size_t n, T* out, const mxx::comm& comm = mxx::comm()) {
    allreduce(in, n, out, std::plus<T>(), comm);
}

template <typename T, typename Func>
inline std::vector<T> allreduce(const T* in, size_t n, Func func, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result(n);
    allreduce(in, n, &result[0], func, comm);
    return result;
}

template <typename T>
inline std::vector<T> allreduce(const T* in, size_t n, const mxx::comm& comm = mxx::comm()) {
    return allreduce(in, n, std::plus<T>(), comm);
}

template <typename T, typename Func>
inline std::vector<T> allreduce(const std::vector<T>& in, Func func, const mxx::comm& comm = mxx::comm()) {
    return allreduce(&in[0], in.size(), func, comm);
}

template <typename T>
inline std::vector<T> allreduce(const std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    return allreduce(in, std::plus<T>(), comm);
}

template <typename T, typename Func>
inline T allreduce(const T& x, Func func, const mxx::comm& comm = mxx::comm()) {
    // get custom op (and type for custom op)
    mxx::custom_op<T> op(std::forward<Func>(func));
    // perform reduction
    T result;
    MPI_Allreduce(const_cast<T*>(&x), &result, 1, op.get_type(), op.get_op(), comm);
    return result;
}

template <typename T>
inline T allreduce(const T& x, const mxx::comm& comm = mxx::comm()) {
    return allreduce(x, std::plus<T>(), comm);
}

/************************
 *  Boolean reductions  *
 ************************/
// useful for testing global conditions, such as termination conditions

inline bool all_of(bool x, const mxx::comm& comm = mxx::comm()) {
    int i = x ? 1 : 0;
    int result;
    MPI_Allreduce(&i, &result, 1, MPI_INT, MPI_LAND, comm);
    return result != 0;
}


inline bool any_of(bool x, const mxx::comm& comm = mxx::comm()) {
    int i = x ? 1 : 0;
    int result;
    MPI_Allreduce(&i, &result, 1, MPI_INT, MPI_LOR, comm);
    return result != 0;
}

inline bool none_of(bool x, const mxx::comm& comm = mxx::comm()) {
    int i = x ? 1 : 0;
    int result;
    MPI_Allreduce(&i, &result, 1, MPI_INT, MPI_LAND, comm);
    return result == 0;
}


template <typename T>
inline bool all_same(const T& x, const mxx::comm& comm = mxx::comm()) {
    T y = mxx::right_shift(x, comm);
    bool same = comm.rank() == 0 || y == x;
    return all_of(same, comm);
}

/*********************************************************************
 *                    Local and Global reductions                    *
 *********************************************************************/

// local reduce

template <typename Iterator, typename Func>
inline typename std::iterator_traits<Iterator>::value_type local_reduce(Iterator begin, Iterator end, Func func) {
    assert(std::distance(begin, end) >= 1);
    typedef typename std::iterator_traits<Iterator>::value_type T;
    T init = std::accumulate(begin+1, end, *begin, func);
    return init;
}

template <typename Iterator>
inline typename std::iterator_traits<Iterator>::value_type local_reduce(Iterator begin, Iterator end) {
    return local_reduce(begin, end, std::plus<typename std::iterator_traits<Iterator>::value_type>());
}

// overloads for std::vector

template <typename T, typename Func>
inline T local_reduce(const std::vector<T>& in, Func func) {
    assert(in.size() >= 1);
    T init = std::accumulate(in.begin()+1, in.end(), in.front(), func);
    return init;
}
template <typename T>
inline T local_reduce(const std::vector<T>& in) {
    return local_reduce(in, std::plus<T>());
}

// global reduce (= local_reduce + allreduce)

template <typename Iterator, typename Func>
inline typename std::iterator_traits<Iterator>::value_type global_reduce(Iterator begin, Iterator end, Func func, const mxx::comm& comm = mxx::comm()) {
    size_t n = std::distance(begin, end);
    typedef typename std::iterator_traits<Iterator>::value_type T;

    if (mxx::allreduce((int)(n >= 1), [](int x, int y) { return (int)(x && y); }, comm) == 0) {
        // some processors have 0 elements
        mxx::comm nonzero_comm = comm.split(n >= 1);
        if (n == 0 && nonzero_comm.size() == comm.size()) {
            // all processes have zero elements, thus return default value
            return T();
        }
        // otherwise reduce only over nonzero processors (subcommunicator)
        int bcast_rank = -1;
        T result;
        if (n >= 1) {
            // local reduction
            result = std::accumulate(begin+1, end, *begin, func);
            // reduction in nonzero subcommunicator
            result = reduce(result, 0, func, nonzero_comm);
            // determine rank of first element of nonzero comm
            if (nonzero_comm.rank() == 0)
                bcast_rank = comm.rank();
        }
        // get rank of processor for bcast
        int bcast_src = mxx::allreduce(bcast_rank, mxx::max<int>(), comm);
        mxx::datatype dt = mxx::get_datatype<T>();
        MPI_Bcast(&result, 1, dt.type(), bcast_src, comm);
        return result;
    } else {
        assert(n >= 1);
        T init = std::accumulate(begin+1, end, *begin, func);
        return allreduce(init, func, comm);
    }
}

template <typename Iterator>
inline typename std::iterator_traits<Iterator>::value_type global_reduce(Iterator begin, Iterator end, const mxx::comm& comm = mxx::comm()) {
    return global_reduce(begin, end, std::plus<typename std::iterator_traits<Iterator>::value_type>(), comm);
}

// overloads for std::vector

template <typename T, typename Func>
inline T global_reduce(const std::vector<T>& in, Func func, const mxx::comm& comm = mxx::comm()) {
    return global_reduce(in.begin(), in.end(), func, comm);
}

template <typename T>
inline T global_reduce(const std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    return global_reduce(in.begin(), in.end(), std::plus<T>(), comm);
}

/*********************************************************************
 *                       max/min with location                       *
 *********************************************************************/

template <typename T>
inline std::pair<T, int> max_element(const T& x, const mxx::comm& comm = mxx::comm()) {
    if (mxx::is_builtin_pair_type<T>::value) {
        MPI_Datatype dt = mxx::datatype_pair<T>::get_type();
        struct {
            T value;
            int rank;
        } in, out;
        in.value = x;
        in.rank = comm.rank();
        MPI_Allreduce(&in, &out, 1, dt, MPI_MAXLOC, comm);
        return std::make_pair(out.value, out.rank);
    } else {
        // use custom operator
        std::pair<T, int> in = std::make_pair(x, comm.rank());
        return mxx::allreduce(in, [](const std::pair<T, int>& x, const std::pair<T, int>& y) { return x.first < y.first ? y : x;}, comm);
    }
}

// vector operation
template <typename T>
inline std::vector<std::pair<T, int>> max_element(const std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    std::vector<std::pair<T, int> > pairin(in.size());
    MXX_ASSERT(mxx::all_same(in.size(), comm));
    for (size_t i = 0; i < in.size(); ++i) {
        pairin[i] = std::make_pair(in[i], comm.rank());
    }
    // don't use MPI_MAXLOC, because it requires re-packing of the data into structs
    return mxx::allreduce(pairin, [](const std::pair<T, int>& x, const std::pair<T, int>& y) { return x.first < y.first ? y : x;}, comm);
}

template <typename T>
inline std::pair<T, int> min_element(const T& x, const mxx::comm& comm = mxx::comm()) {
    if (mxx::is_builtin_pair_type<T>::value) {
        MPI_Datatype dt = mxx::datatype_pair<T>::get_type();
        struct {
            T value;
            int rank;
        } in, out;
        in.value = x;
        in.rank = comm.rank();
        MPI_Allreduce(&in, &out, 1, dt, MPI_MINLOC, comm);
        return std::make_pair(out.value, out.rank);
    } else {
        // use custom operator
        std::pair<T, int> in = std::make_pair(x, comm.rank());
        return mxx::allreduce(in, [](const std::pair<T, int>& x, const std::pair<T, int>& y) { return x.first > y.first ? y : x;}, comm);
    }
}

// vector operation
template <typename T>
inline std::vector<std::pair<T, int>> min_element(const std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    std::vector<std::pair<T, int> > pairin(in.size());
    MXX_ASSERT(mxx::all_same(in.size(), comm));
    for (size_t i = 0; i < in.size(); ++i) {
        pairin[i] = std::make_pair(in[i], comm.rank());
    }
    // don't use MPI_MAXLOC, because it requires re-packing of the data into structs
    return mxx::allreduce(pairin, [](const std::pair<T, int>& x, const std::pair<T, int>& y) { return x.first > y.first ? y : x;}, comm);
}


/*********************************************************************
 *                               Scan                                *
 *********************************************************************/

// reduce over vectors

template <typename T, typename Func>
inline void scan_vec(const T* in, size_t n, T* out, Func func, const mxx::comm& comm = mxx::comm()) {
    // get op
    mxx::custom_op<T> op(std::forward<Func>(func));
    MPI_Scan(const_cast<T*>(in), out, n, op.get_type(), op.get_op(), comm);
}
template <typename T, typename Func>
inline std::vector<T> scan_vec(const T* in, size_t n, Func func, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result(n);
    scan_vec(in, n, &result[0], func, comm);
    return result;
}
template <typename T, typename Func>
inline std::vector<T> scan_vec(const std::vector<T>& x, Func func, const mxx::comm& comm = mxx::comm()) {
    return scan_vec(&x[0], x.size(), func, comm);
}

// single element per processor

template <typename T, typename Func>
inline T scan(const T& x, Func func, const mxx::comm& comm = mxx::comm()) {
    // get op
    mxx::custom_op<T> op(std::forward<Func>(func));
    T result;
    MPI_Scan(const_cast<T*>(&x), &result, 1, op.get_type(), op.get_op(), comm);
    return result;
}

template <typename T>
inline T scan(const T& x, const mxx::comm& comm = mxx::comm()) {
    return scan(x, std::plus<T>(), comm);
}

// local scan

template <typename InIterator, typename OutIterator, typename Func>
void local_scan(InIterator begin, InIterator end, OutIterator out, Func func) {
    // return if there's nothing here
    if (std::distance(begin, end) == 0)
        return;
    typedef typename std::iterator_traits<OutIterator>::value_type T;

    // start from first element
    T val = *begin;
    *out = *begin;
    ++begin;
    ++out;

    // calculate the inclusive prefix sum
    while (begin != end) {
        val = func(val,*begin);
        *out = val;
        ++begin;
        ++out;
    }
}


// inplace!
template <typename Iterator, typename Func>
void local_scan_inplace(Iterator begin, Iterator end, Func func) {
    // return if there's nothing here
    if (std::distance(begin, end) == 0)
        return;
    typedef typename std::iterator_traits<Iterator>::value_type T;

    // start from first element
    T val = *begin;
    ++begin;

    // calculate the inclusive prefix sum
    while (begin != end) {
        val = func(val,*begin);
        *begin = val;
        ++begin;
    }
}

template <typename InIterator, typename OutIterator>
inline void local_scan(InIterator begin, InIterator end, OutIterator out) {
    return local_scan(begin, end, out, std::plus<typename std::iterator_traits<OutIterator>::value_type>());
}

template <typename Iterator>
inline void local_scan_inplace(Iterator begin, Iterator end) {
    return local_scan_inplace(begin, end, std::plus<typename std::iterator_traits<Iterator>::value_type>());
}

// std::vector overloads
template <typename T, typename Func>
inline void local_scan_inplace(std::vector<T>& in, Func func) {
    local_scan_inplace(in.begin(), in.end(), func);
}

template <typename T>
inline void local_scan_inplace(std::vector<T>& in) {
    local_scan_inplace(in.begin(), in.end(), std::plus<T>());
}

template <typename T, typename Func>
inline std::vector<T> local_scan(const std::vector<T>& in, Func func) {
    std::vector<T> result(in.size());
    local_scan(in.begin(), in.end(), result.begin(), func);
    return result;
}

template <typename T>
inline std::vector<T> local_scan(const std::vector<T>& in) {
    std::vector<T> result(in.size());
    local_scan(in.begin(), in.end(), result.begin(), std::plus<T>());
    return result;
}

// global scans

template <typename InIterator, typename OutIterator, typename Func>
void global_scan(InIterator begin, InIterator end, OutIterator out, Func func, const mxx::comm& comm = mxx::comm()) {
    OutIterator o = out;
    size_t n = std::distance(begin, end);
    // create subcommunicator for those processes which contain elements
    mxx::comm nonzero_comm = comm.split(n > 0);
    if (n > 0) {
        // local scan
        local_scan(begin, end, out, func);
        // mxx::scan
        typedef typename std::iterator_traits<OutIterator>::value_type T;
        T sum = T();
        if (n > 0)
            sum = *(out+(n-1));
        T presum = exscan(sum, func, nonzero_comm);
        // accumulate previous sum on all local elements
        for (size_t i = 0; i < n; ++i) {
            *o = func(presum, *o);
            ++o;
        }
    }
}


// inplace!
template <typename Iterator, typename Func>
inline void global_scan_inplace(Iterator begin, Iterator end, Func func, const mxx::comm& comm = mxx::comm()) {
    Iterator o = begin;
    size_t n = std::distance(begin, end);
    mxx::comm nonzero_comm = comm.split(n > 0);
    if (n > 0) {
        // local inplace scan
        local_scan_inplace(begin, end, func);
        // mxx::exscan
        typedef typename std::iterator_traits<Iterator>::value_type T;
        T sum = *(begin + (n-1));
        T presum = exscan(sum, func, nonzero_comm);

        // accumulate previous sum on all local elements
        for (size_t i = 0; i < n; ++i) {
            *o = func(presum, *o);
            ++o;
        }
    }
}

template <typename InIterator, typename OutIterator>
inline void global_scan(InIterator begin, InIterator end, OutIterator out, const mxx::comm& comm = mxx::comm()) {
    return global_scan(begin, end, out, std::plus<typename std::iterator_traits<OutIterator>::value_type>(), comm);
}

template <typename Iterator>
inline void global_scan_inplace(Iterator begin, Iterator end, const mxx::comm& comm = mxx::comm()) {
    return global_scan_inplace(begin, end, std::plus<typename std::iterator_traits<Iterator>::value_type>(), comm);
}

// std::vector overloads
template <typename T, typename Func>
inline void global_scan_inplace(std::vector<T>& in, Func func, const mxx::comm& comm = mxx::comm()) {
    global_scan_inplace(in.begin(), in.end(), func, comm);
}

template <typename T>
inline void global_scan_inplace(std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    global_scan_inplace(in.begin(), in.end(), std::plus<T>(), comm);
}

template <typename T, typename Func>
inline std::vector<T> global_scan(const std::vector<T>& in, Func func, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result(in.size());
    global_scan(in.begin(), in.end(), result.begin(), func, comm);
    return result;
}

template <typename T>
inline std::vector<T> global_scan(const std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result(in.size());
    global_scan(in.begin(), in.end(), result.begin(), std::plus<T>(), comm);
    return result;
}

/*********************************************************************
 *                              Exscan                               *
 *********************************************************************/

// reduce over vectors

template <typename T, typename Func>
inline void exscan_vec(const T* in, size_t n, T* out, Func func, const mxx::comm& comm = mxx::comm()) {
    // get op
    mxx::custom_op<T> op(std::forward<Func>(func));
    MPI_Exscan(const_cast<T*>(in), out, n, op.get_type(), op.get_op(), comm);
}
template <typename T, typename Func>
inline std::vector<T> exscan_vec(const T* in, size_t n, Func func, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result(n, T());
    exscan_vec(in, n, &result[0], func, comm);
    return result;
}
template <typename T, typename Func>
inline std::vector<T> exscan_vec(const std::vector<T>& x, Func func, const mxx::comm& comm = mxx::comm()) {
    return exscan_vec(&x[0], x.size(), func, comm);
}


// single element

template <typename T, typename Func>
T exscan(const T& x, Func func, const mxx::comm& comm = mxx::comm()) {
    // get op
    mxx::custom_op<T> op(std::forward<Func>(func));
    // perform reduction
    T result;
    MPI_Exscan(const_cast<T*>(&x), &result, 1, op.get_type(), op.get_op(), comm);
    if (comm.rank() == 0)
      result = T();
    return result;
}

template <typename T>
T exscan(const T& x, const mxx::comm& comm = mxx::comm()) {
    return exscan(x, std::plus<T>(), comm);
}


// local exscan

template <typename InIterator, typename OutIterator, typename Func>
void local_exscan(InIterator begin, InIterator end, OutIterator out, Func func) {
    // return if there's nothing here
    if (std::distance(begin, end) == 0)
        return;
    typedef typename std::iterator_traits<OutIterator>::value_type T;

    // start from first element
    T val = *begin;
    *out = T();
    ++begin;
    ++out;

    // calculate the exclusive prefix sum
    while (begin != end) {
        T tmp = val;
        val = func(val,*begin);
        *out = tmp;
        ++begin;
        ++out;
    }
}


// inplace!
template <typename Iterator, typename Func>
void local_exscan_inplace(Iterator begin, Iterator end, Func func) {
    // return if there's nothing here
    if (std::distance(begin, end) == 0)
        return;
    typedef typename std::iterator_traits<Iterator>::value_type T;

    // start from first element
    T val = *begin;
    *begin = T();
    ++begin;

    // calculate the exclusive prefix sum
    while (begin != end) {
        T tmp = val;
        val = func(val,*begin);
        *begin = tmp;
        ++begin;
    }
}

template <typename InIterator, typename OutIterator>
inline void local_exscan(InIterator begin, InIterator end, OutIterator out) {
    return local_exscan(begin, end, out, std::plus<typename std::iterator_traits<OutIterator>::value_type>());
}

template <typename Iterator>
inline void local_exscan_inplace(Iterator begin, Iterator end) {
    return local_exscan_inplace(begin, end, std::plus<typename std::iterator_traits<Iterator>::value_type>());
}

// std::vector overloads
template <typename T, typename Func>
inline void local_exscan_inplace(std::vector<T>& in, Func func) {
    local_exscan_inplace(in.begin(), in.end(), func);
}

template <typename T>
inline void local_exscan_inplace(std::vector<T>& in) {
    local_exscan_inplace(in.begin(), in.end(), std::plus<T>());
}

template <typename T, typename Func>
inline std::vector<T> local_exscan(const std::vector<T>& in, Func func) {
    std::vector<T> result(in.size());
    local_exscan(in.begin(), in.end(), result.begin(), func);
    return result;
}

template <typename T>
inline std::vector<T> local_exscan(const std::vector<T>& in) {
    std::vector<T> result(in.size());
    local_exscan(in.begin(), in.end(), result.begin(), std::plus<T>());
    return result;
}

// global scans

template <typename InIterator, typename OutIterator, typename Func>
void global_exscan(InIterator begin, InIterator end, OutIterator out, Func func, const mxx::comm& comm = mxx::comm()) {
    OutIterator o = out;
    size_t n = std::distance(begin, end);
    mxx::comm nonzero_comm = comm.split(n > 0);
    if (n > 0) {
        typedef typename std::iterator_traits<OutIterator>::value_type T;
        T sum = *(begin+(n-1));
        // local scan
        local_exscan(begin, end, out, func);
        // mxx::scan
        sum += *(o + (n-1));
        T presum = exscan(sum, func, nonzero_comm);
        *o++ = presum;
        // accumulate previous sum on all local elements
        for (size_t i = 1; i < n; ++i) {
            *o = func(presum, *o);
            ++o;
        }
    }
}


// inplace!
template <typename Iterator, typename Func>
inline void global_exscan_inplace(Iterator begin, Iterator end, Func func, const mxx::comm& comm = mxx::comm()) {
    Iterator o = begin;
    size_t n = std::distance(begin, end);
    mxx::comm nonzero_comm = comm.split(n > 0);
    if (n > 0) {
        typedef typename std::iterator_traits<Iterator>::value_type T;
        T sum = *(begin+(n-1));
        // local inplace scan
        local_exscan_inplace(begin, end, func);
        sum += *(begin+(n-1));
        // mxx::exscan
        T presum = exscan(sum, func, nonzero_comm);
        *o++ = presum;
        // accumulate previous sum on all local elements
        for (size_t i = 1; i < n; ++i) {
            *o = func(presum, *o);
            ++o;
        }
    }
}

template <typename InIterator, typename OutIterator>
inline void global_exscan(InIterator begin, InIterator end, OutIterator out, const mxx::comm& comm = mxx::comm()) {
    return global_exscan(begin, end, out, std::plus<typename std::iterator_traits<OutIterator>::value_type>(), comm);
}

template <typename Iterator>
inline void global_exscan_inplace(Iterator begin, Iterator end, const mxx::comm& comm = mxx::comm()) {
    return global_exscan_inplace(begin, end, std::plus<typename std::iterator_traits<Iterator>::value_type>(), comm);
}

// std::vector overloads
template <typename T, typename Func>
inline void global_exscan_inplace(std::vector<T>& in, Func func, const mxx::comm& comm = mxx::comm()) {
    global_exscan_inplace(in.begin(), in.end(), func, comm);
}

template <typename T>
inline void global_exscan_inplace(std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    global_exscan_inplace(in.begin(), in.end(), std::plus<T>(), comm);
}

template <typename T, typename Func>
inline std::vector<T> global_exscan(const std::vector<T>& in, Func func, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result(in.size());
    global_exscan(in.begin(), in.end(), result.begin(), func, comm);
    return result;
}

template <typename T>
inline std::vector<T> global_exscan(const std::vector<T>& in, const mxx::comm& comm = mxx::comm()) {
    std::vector<T> result(in.size());
    global_exscan(in.begin(), in.end(), result.begin(), std::plus<T>(), comm);
    return result;
}

/****************************************************
 *  reverse reductions (with reverse communicator)  *
 ****************************************************/

/*
inline void rev_comm(MPI_Comm comm, MPI_Comm& rev)
{
    // get MPI parameters
    int rank;
    int p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);
    MPI_Comm_split(comm, 0, p - rank, &rev);
}

template <typename T>
T reverse_exscan(const T& x, const mxx::comm& comm = mxx::comm()) {
    MPI_Comm rev;
    rev_comm(comm, rev);
    T result = exscan(x, rev);
    MPI_Comm_free(&rev);
    return result;
}

template <typename T, typename Func>
T reverse_exscan(const T& x, Func func, const mxx::comm& comm = mxx::comm()) {
    MPI_Comm rev;
    rev_comm(comm, rev);
    T result = exscan(x, func, rev);
    MPI_Comm_free(&rev);
    return result;
}

template <typename T>
T reverse_scan(const T& x, const mxx::comm& comm = mxx::comm()) {
    MPI_Comm rev;
    rev_comm(comm, rev);
    T result = scan(x, rev);
    MPI_Comm_free(&rev);
    return result;
}

template <typename T, typename Func>
T reverse_scan(const T& x, Func func, const mxx::comm& comm = mxx::comm()) {
    MPI_Comm rev;
    rev_comm(comm, rev);
    T result = scan(x, func, rev);
    MPI_Comm_free(&rev);
    return result;
}
*/


/*********************
 *  Specialized ops  *
 *********************/


} // namespace mxx

#define MXX_REDUCTION_DONE
#include "comm_def.hpp"

#endif // MXX_REDUCTION_HPP
