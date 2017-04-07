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
 * @file    datatypes.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   MPI Datatypes for C++ types.
 */

#ifndef MXX_DATATYPES_HPP
#define MXX_DATATYPES_HPP

// MPI include
#include <mpi.h>

// C++ includes
#include <vector>
#include <map>
#include <array>
#include <tuple>
#include <numeric>
#include <limits>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <iostream>

#include "common.hpp"
#include "type_traits.hpp"



namespace mxx
{

/*
 * Mapping of C/C++ types to MPI datatypes.
 *
 * Possible ways of implementation
 * 1) templated function get_mpi_dt<template T>(); (e.g. boost)
 *     -> doesn't allow partial template specialization
 *        (needed for std::pair, std::tuple etc)
 * 2) overloaded functions with type deduction
 *     -> can't properly clean up types (MPI_Type_free)
 * 3) templated class with static method (-> allows partial specialization)
 * 4) templated class with member method (allows proper C++ like Type_free)!!!
 *
 * Using option 4 due to best fit.
 *
 * TODO: other (maybe better) possibility:
 * 5) using templated class with static methods and cache all created datatypes
 *    in global static map (typeid(T) -> MPI_Datatype), freeing upon global
 *    destruction
 */

// TODO:
// - [x] compile time checker for builtin types
// - [ ] put attr_map in here!
// - [ ] add to-string and caching
// - [ ] implement MPI type introspection (get envelope)
// - [ ] (see wrappers in official MPI C++ bindings)

template <typename T>
class is_builtin_type : public std::false_type {};

class datatype;

// TODO: there has to be a better way for this
} // namespace mxx

std::false_type make_datatype();

namespace mxx {

template <typename T>
class static_datatype_builder;

template <typename T>
datatype get_datatype();

/*
template <typename T>
datatype build_datatype(const T&);
*/

template <typename T, typename... Args>
datatype built_custom_datatype(T*, Args&...);

// MPI_Datatype wrapper
class datatype {
protected:
    /*
    template <typename T>
    friend datatype get_datatype();
    */

    template <typename T>
    friend datatype build_datatype(const T&);

    template <typename T, typename... Args>
    friend datatype built_custom_datatype(T*, Args&...);
    template <typename T>
    friend class custom_datatype_builder;
    template <typename T, typename Derived>
    friend class datatype_builder_base;
public:
    datatype() : mpitype(MPI_DATATYPE_NULL), builtin(true) {
    }

    // TODO: make this constructor protected and add friend access to builder function
    datatype(MPI_Datatype mpidt, bool builtin) : mpitype(mpidt), builtin(builtin) {
    }

    // copy constructor
    datatype(const datatype& o) = delete;

    // move constructor
    datatype(datatype&& o) {
        builtin = o.builtin;
        mpitype = o.mpitype;
        o.mpitype = MPI_DATATYPE_NULL;
        o.builtin = true;
    }

    // copy assignment
    datatype& operator=(const datatype& o) = delete;

    // move assignment
    datatype& operator=(datatype&& o) {
        builtin = o.builtin;
        mpitype = o.mpitype;
        o.mpitype = MPI_DATATYPE_NULL;
        o.builtin = true;
        return *this;
    }

    MPI_Datatype type() const {
        return mpitype;
    }

    datatype vector(size_t count, size_t blocklen, size_t stride) const {
        // TODO what if any of the parameters > MAX_INT? -> use Type_create_struct instead!
        //  -> if blocklen >= MAX_INT then use contiguous first, then struct with count elements
        datatype result;
        result.builtin = false;
        MPI_Type_vector(count, blocklen, stride, this->mpitype, &result.mpitype);
        return result;
    }

    datatype contiguous(size_t count) const {
        datatype result;
        if (count <= mxx::max_int) {
            result.builtin = false;
            MPI_Type_contiguous(count, this->mpitype, &result.mpitype);
            MPI_Type_commit(&result.mpitype);
        } else {
            result.builtin = false;
            // create custom data types of blocks and remainder
            std::size_t intmax = mxx::max_int;
            std::size_t nblocks = count / intmax;
            std::size_t rem = count % intmax;

            // create block and remainder data types
            MPI_Datatype _block;
            MPI_Type_contiguous(mxx::max_int, this->mpitype, &_block);
            MPI_Datatype _blocks;
            // create two contiguous types for blocks and remainder
            MPI_Type_contiguous(nblocks, _block, &_blocks);

            if (rem > 0) {
                MPI_Datatype _remainder;
                MPI_Type_contiguous(rem, this->mpitype, &_remainder);

                // create struct for the concatenation of this type
                MPI_Aint lb, extent;
                MPI_Type_get_extent(this->mpitype, &lb, &extent);
                MPI_Aint displ = nblocks*intmax*extent;
                MPI_Aint displs[2] = {0, displ};
                int blocklen[2] = {1, 1};
                MPI_Datatype mpitypes[2] = {_blocks, _remainder};
                MPI_Type_create_struct(2, blocklen, displs, mpitypes, &result.mpitype);
                MPI_Type_commit(&result.mpitype);
                MPI_Type_free(&_remainder);
                MPI_Type_free(&_blocks);
            } else {
                result.mpitype = _blocks;
                MPI_Type_commit(&result.mpitype);
            }
        }
        return result;
    }

    std::pair<MPI_Aint, MPI_Aint> get_extent() {
        std::pair<MPI_Aint, MPI_Aint> e;
        MPI_Type_get_extent(mpitype, &e.first, &e.second);
        return e;
    }

    std::pair<MPI_Aint, MPI_Aint> get_true_extent() {
        std::pair<MPI_Aint, MPI_Aint> e;
        MPI_Type_get_true_extent(mpitype, &e.first, &e.second);
        return e;
    }

    // TODO: get envelope + get_contents for printing of datatype!
    void get_envelope() {
    }

    // TODO: indexed/hindexed?

    virtual ~datatype() {
        // TODO: don't free, but decrease counter in type cache
        if (!builtin)
            MPI_Type_free(&mpitype);
    }
private:
    MPI_Datatype mpitype;
    bool builtin;
};


// defined generalized datatype builder
template <typename T>
struct datatype_builder {};

/*********************************************************************
 *                     Define built-in datatypes                     *
 *********************************************************************/

#define MXX_DATATYPE_MPI_BUILTIN(ctype, mpi_type)                           \
template <> struct datatype_builder<ctype> {                                \
    static MPI_Datatype get_type() {return mpi_type;}                       \
    static size_t num_basic_elements() { return 1;}                         \
};                                                                          \
                                                                            \
template <> class is_builtin_type<ctype> : public std::true_type {};        \


// calls the given macro on each pair of builtin type and corresponding
// MPI_Datatype
#define MXX_FOR_ALL_BUILTIN(BUILTIN_TYPE)                                      \
/* char */                                                                     \
BUILTIN_TYPE(char, MPI_CHAR);                                                  \
BUILTIN_TYPE(unsigned char, MPI_UNSIGNED_CHAR);                                \
BUILTIN_TYPE(signed char, MPI_SIGNED_CHAR);                                    \
                                                                               \
/* short */                                                                    \
BUILTIN_TYPE(unsigned short, MPI_UNSIGNED_SHORT);                              \
BUILTIN_TYPE(short, MPI_SHORT);                                                \
                                                                               \
/* int */                                                                      \
BUILTIN_TYPE(unsigned int, MPI_UNSIGNED);                                      \
BUILTIN_TYPE(int, MPI_INT);                                                    \
                                                                               \
/* long */                                                                     \
BUILTIN_TYPE(unsigned long, MPI_UNSIGNED_LONG);                                \
BUILTIN_TYPE(long, MPI_LONG);                                                  \
                                                                               \
/* long long */                                                                \
BUILTIN_TYPE(unsigned long long, MPI_UNSIGNED_LONG_LONG);                      \
BUILTIN_TYPE(long long, MPI_LONG_LONG);                                        \
                                                                               \
/* floats */                                                                   \
BUILTIN_TYPE(float, MPI_FLOAT);                                                \
BUILTIN_TYPE(double, MPI_DOUBLE);                                              \
BUILTIN_TYPE(long double, MPI_LONG_DOUBLE);                                    \


MXX_FOR_ALL_BUILTIN(MXX_DATATYPE_MPI_BUILTIN);

#undef MXX_DATATYPE_MPI_BUILTIN


struct datatype_name {
    std::string mpi_name;
    std::string c_name;
    std::string typeid_name;

    datatype_name(const std::string& mpi_name, const std::string& c_name, const std::string& typeid_name)
      : mpi_name(mpi_name), c_name(c_name), typeid_name(typeid_name) {}

    datatype_name() {}
    datatype_name(const datatype_name& o) = default;
    datatype_name(datatype_name&& o) = default;
};

inline std::ostream& operator<<(std::ostream& os, const datatype_name& n) {
    return os << "(" << n.mpi_name << "," << n.c_name << "," << n.typeid_name << ")";
}

// define reverse mapping of datatypes for type decoding
#define MXX_INSERT_NAME_INTO_MAP(ctype, mpi_type) \
m.emplace(mpi_type, datatype_name(#mpi_type, #ctype, typeid(ctype).name()))

class builtin_typename_map {
private:
    static std::unordered_map<MPI_Datatype, datatype_name> init_map() {
        std::unordered_map<MPI_Datatype, datatype_name> m;
        MXX_FOR_ALL_BUILTIN(MXX_INSERT_NAME_INTO_MAP);
        return m;
    }

public:
    static std::unordered_map<MPI_Datatype, datatype_name>& get_map() {
        // C++11 standard guarantuess that a static variable gets instantiated
        // in a threadsafe manner
        static std::unordered_map<MPI_Datatype, datatype_name> m = init_map();
        return m;
    }
    static std::string get_typeid_name(const MPI_Datatype& dt) {
        return get_map()[dt].typeid_name;
    }
    static std::string get_c_name(const MPI_Datatype& dt) {
        return get_map()[dt].c_name;
    }
    static std::string get_mpi_name(const MPI_Datatype& dt) {
      return get_map()[dt].mpi_name;
    }
};



/*********************************************************************
 *                 Pair types for MINLOC and MAXLOC                  *
 *********************************************************************/

template <typename T>
struct datatype_pair {
    static MPI_Datatype get_type() {
        return MPI_DATATYPE_NULL;
    }
};

template <typename T>
class is_builtin_pair_type : public std::false_type {};

#define MXX_DATATYPE_BUILTIN_PAIR(ctype, mpi_type)                          \
template <> struct datatype_pair<ctype> {                                   \
    static MPI_Datatype get_type() {                                        \
        return mpi_type;                                                    \
    }                                                                       \
};                                                                          \
template <> class is_builtin_pair_type<ctype> : public std::true_type {};   \

// integers-integer pairs
MXX_DATATYPE_BUILTIN_PAIR(short, MPI_SHORT_INT);
MXX_DATATYPE_BUILTIN_PAIR(int, MPI_2INT);
MXX_DATATYPE_BUILTIN_PAIR(long, MPI_LONG_INT);

// floats
MXX_DATATYPE_BUILTIN_PAIR(float, MPI_FLOAT_INT);
MXX_DATATYPE_BUILTIN_PAIR(double, MPI_DOUBLE_INT);
MXX_DATATYPE_BUILTIN_PAIR(long double, MPI_LONG_DOUBLE_INT);


#undef MXX_DATATYPE_BUILTIN_PAIR

template <typename T>
struct has_datatype<T, typename std::enable_if<mxx::is_builtin_type<T>::value>::type> : std::true_type {};

// TODO: extend this!
template <typename T, typename U>
struct has_datatype<std::pair<T, U>> : std::true_type {};


/**
 * @brief   MPI datatype mapping for std::array
 */
template <typename T, std::size_t size>
struct datatype_builder<std::array<T, size> > {
    static MPI_Datatype get_type() {
        MPI_Datatype _type;
        MPI_Datatype base_type = get_datatype<T>().type();
        MPI_Type_contiguous(size, base_type, &_type);
        MPI_Type_commit(&_type);
        return _type;
    }
    static size_t num_basic_elements() {
        return size*datatype_builder<T>::num_basic_elements();
    }
};

/**
 * @brief   MPI datatype mapping for std::pair
 */
template <typename T1, typename T2>
struct datatype_builder<std::pair<T1, T2> > {
    static MPI_Datatype get_type() {
        MPI_Datatype _type;

        int blocklen[2] = {1, 1};
        MPI_Aint displs[2] = {0,0};
        // get actual displacement (in case of padding in the structure)
        std::pair<T1, T2> p;
        MPI_Aint p_adr, t1_adr, t2_adr;
        MPI_Get_address(&p, &p_adr);
        MPI_Get_address(&p.first, &t1_adr);
        MPI_Get_address(&p.second, &t2_adr);
        displs[0] = t1_adr - p_adr;
        displs[1] = t2_adr - p_adr;

        // create type
        // TODO: use cached type!
        MPI_Datatype types[2] = {datatype_builder<T1>::get_type(), datatype_builder<T2>::get_type()};
        // in case elements are represented the opposite way around in
        // the pair (gcc does so), then swap them
        if (displs[0] > displs[1]) {
            std::swap(displs[0], displs[1]);
            std::swap(types[0], types[1]);
        }
        // create MPI_Datatype (resized to actual sizeof())
        MPI_Datatype struct_type;
        MPI_Type_create_struct(2, blocklen, displs, types, &struct_type);
        MPI_Type_create_resized(struct_type, 0, sizeof(p), &_type);
        MPI_Type_commit(&_type);
        MPI_Type_free(&struct_type);
        return _type;
    }

    static size_t num_basic_elements() {
        return datatype_builder<T1>::num_basic_elements() + datatype_builder<T2>::num_basic_elements();
    }
};


// fill in MPI types
template <std::size_t N, std::size_t I, class ...Types>
struct tuple_members {
    static void get(std::map<MPI_Aint, MPI_Datatype>& members) {
        // init tuple to get measurement offsets
        // TODO: use null-ref instead of actual instantiation
        std::tuple<Types...> tuple;

        // get member displacement
        MPI_Aint t_adr, elem_adr;
        MPI_Get_address(&tuple, &t_adr);
        MPI_Get_address(&std::get<N-I>(tuple), &elem_adr);
        // byte offset from beginning of tuple
        MPI_Aint displ = elem_adr - t_adr;
        // fill in type
        // TODO: use cached type!?
        MPI_Datatype mpi_dt = datatype_builder<typename std::tuple_element<N-I,std::tuple<Types...>>::type>::get_type(); //std::get<N-I>(datatypes).type();

        // add to map
        members[displ] = mpi_dt;

        // recursively (during compile time) call same function
        tuple_members<N,I-1, Types...>::get(members);
    }
};

// Base case of meta-recursion
template <std::size_t N, class ...Types>
struct tuple_members<N, 0, Types...> {
    static void get(std::map<MPI_Aint, MPI_Datatype>&) {
    }
};

template <class...Types>
struct tuple_basic_els;

template <class T, class...Types>
struct tuple_basic_els<T,Types...>
{
    static size_t get_num() {
        return datatype_builder<T>::num_basic_elements() + tuple_basic_els<Types...>::get_num();
    }
};

template <class T>
struct tuple_basic_els<T>
{
    static size_t get_num() {
        return datatype_builder<T>::num_basic_elements();
    }
};

/**
 * @brief   MPI datatype mapping for std::tuple
 */
template <class ...Types>
struct datatype_builder<std::tuple<Types...> > {
  // tuple type
  typedef std::tuple<Types...> tuple_t;
  // number of elements of the typle
  static constexpr std::size_t size = std::tuple_size<tuple_t>::value;

  /// returns the MPI_Datatype for the tuple
  static MPI_Datatype get_type() {
        MPI_Datatype _type;
        // fill in the block lengths to 1 each
        int blocklen[size];
        for (std::size_t i = 0; i < size; ++i) {
            blocklen[i] = 1;
        }

        // get the member displacement and type info for the tuple using
        // meta-recursion
        std::map<MPI_Aint, MPI_Datatype> members;
        tuple_members<size,size,Types...>::get(members);

        // fill displacements and types according to in-memory order in tuple
        // NOTE: the in-memory order is not necessarily the same as the order
        // of types as accessed by std::get
        // For gcc the order is actually reversed!
        // Hence, we use a std::map to collect the order information prior
        // to creating the displacement and type arrays
        std::array<MPI_Aint, size> displs;
        std::array<MPI_Datatype, size> mpitypes;
        std::size_t i = 0;
        for (std::map<MPI_Aint, MPI_Datatype>::iterator it = members.begin();
             it != members.end(); ++it) {
            displs[i] = it->first;
            mpitypes[i] = it->second;
            ++i;
        }

        // create type
        MPI_Datatype struct_type;
        MPI_Type_create_struct(size, blocklen, &displs[0], &mpitypes[0], &struct_type);
        MPI_Type_create_resized(struct_type, 0, sizeof(tuple_t), &_type);
        MPI_Type_commit(&_type);
        MPI_Type_free(&struct_type);

        return _type;
    }

    static size_t num_basic_elements() {
        return tuple_basic_els<Types...>::get_num();
    }
};

template <typename Derived>
class recursive_processor {
private:
    template <typename M>
    void process_one(M&& m) {
        static_cast<Derived*>(this)->process(std::forward<M>(m));
    }

public:
    // end of recursion
    void process() {}

    template <typename M, typename... Members>
    void process(M&& m, Members&&...vargs) {
        process_one(std::forward<M>(m));
        process(std::forward<Members>(vargs)...);
    }

    // the call operator processes everything
    template <typename... Members>
    void operator()(Members&&...vargs) {
      process(std::forward<Members>(vargs)...);
    }
};

template <typename T, typename Derived>
class datatype_builder_base {
    // saves information about the members (displacement + MPI_Datatype)
    std::map<MPI_Aint, ::mxx::datatype> members;
public:
    template <typename M>
    void add_member_by_offset(size_t offset) {
        // get the underlying datatype
        datatype dt = ::mxx::get_datatype<M>();

        MPI_Aint displ = offset;
        // assert this is actually a member of the type T
        MXX_ASSERT(0 <= displ && displ + sizeof(M) <= sizeof(T));
        // add to map
        members[displ] = std::move(dt);
    }

    // returns the datatype for all added members
    datatype get_datatype() const {
        MXX_ASSERT(members.size() > 0);
        // create the blocklength, displacements, and datatype arrays
        size_t n_members = members.size();
        std::vector<MPI_Aint> displs(n_members);
        std::vector<MPI_Datatype> mpitypes(n_members);
        std::vector<int> blen(n_members, 1);
        std::size_t i = 0;
        for (std::map<MPI_Aint, ::mxx::datatype>::const_iterator it = members.begin();
                it != members.end(); ++it) {
            displs[i] = it->first;
            mpitypes[i] = it->second.type();
            ++i;
        }

        // create type
        MPI_Datatype _type;
        MPI_Datatype struct_type;
        MPI_Type_create_struct(n_members, &blen[0], &displs[0], &mpitypes[0], &struct_type);
        MPI_Type_create_resized(struct_type, 0, sizeof(T), &_type);
        MPI_Type_commit(&_type);
        MPI_Type_free(&struct_type);

        // return mxx::datatype wrapper
        datatype dt(_type, false);
        return dt;
    }
};

template <typename T>
class value_datatype_builder : public datatype_builder_base<T, value_datatype_builder<T>>, recursive_processor<value_datatype_builder<T>> {
private:
    // reference to the type we're building the custom datatype for
    const T& that;
    typedef datatype_builder_base<T, value_datatype_builder<T>> base_type;
public:

    value_datatype_builder(const T& value) : base_type(), that(value) {}

    // custom add_member function which adds members by their offset to `&that`
    template <typename M>
    void add_member(const M& member) {
        // get member displacement
        MPI_Aint t_adr, elem_adr;

        MPI_Get_address((void*)&that, &t_adr);
        MPI_Get_address((void*)&member, &elem_adr);

        // byte offset from beginning of tuple
        MPI_Aint displ = elem_adr - t_adr;
        this->template add_member_by_offset<M>(displ);
    }

    template <typename M>
    void process(M&& m) {
        add_member(std::forward<T>(m));
    }
};

// determine the offset of a `pointer to member` type without instantiation
template <typename T, typename Base, typename M>
typename std::enable_if<std::is_base_of<Base, T>::value, size_t>::type
offset_of(M Base::* m) {
    return reinterpret_cast<size_t>(&(((T*)nullptr)->*m));
}

template <typename T>
class static_datatype_builder : public datatype_builder_base<T, static_datatype_builder<T>> {
private:
    typedef datatype_builder_base<T, static_datatype_builder<T>> base_type;
public:
    // add members via "pointer to member" types
    template <typename M>
    void add_member(M T::*m) {
        this->template add_member_by_offset<M>(offset_of<T, T, M>(m));
    }

    // support adding members of base classes
    template <typename Base, typename M>
    typename std::enable_if<std::is_base_of<Base, T>::value, void>::type
    add_member(M Base::*m) {
        this->template add_member_by_offset<M>(offset_of<T, Base, M>(m));
    }
};


/*
 * "templates" for different kinds of data structures.
 * Inherit from these to specialize for your own type easily.
 */


/**
 * @brief   A contiguous datatype of the same base type
 */
template <typename T, std::size_t size>
struct datatype_contiguous {
    static_assert(size <= std::numeric_limits<int>::max(),
                  "Compile time contiguous types only support sizes up to INT_MAX");
    static MPI_Datatype get_type() {
        MPI_Datatype _type;
        datatype _base_type = get_datatype<T>();
        MPI_Type_contiguous(size, _base_type.type(), &_type);
        MPI_Type_commit(&_type);
        return _type;
    }
    static size_t num_basic_elements() {
        return size*datatype_builder<T>::num_basic_elements();
    }
};

/*
 * Runtime selection of size
 */
/*
template <typename T>
struct datatype_contiguous<T,0> {
    static MPI_Datatype get_type(size_t size) {
        datatype dt = get_datatype<T>().contiguous(size);
        MPI_Datatype mpidt;
        MPI_Type_dup(dt.type(), &mpidt);
        MPI_Type_commit(&mpidt);
        return mpidt;
    }
};
*/

MXX_DEFINE_IS_GLOBAL_FUNC(make_datatype)
MXX_DEFINE_HAS_STATIC_MEMBER(get_type)

template <typename T>
struct has_builder : has_static_member_get_type<datatype_builder<T>, MPI_Datatype()> {};

// basically <=> `has_datatype` (TODO consistent naming)
template <typename T, typename Enable = void>
struct is_trivial_type : std::false_type {};

template <typename T>
struct is_trivial_type<T, typename std::enable_if<
 has_static_member_datatype<T, void(static_datatype_builder<T>&)>::value
 || has_member_datatype<T, void(value_datatype_builder<T>&)>::value
 || is_global_func_make_datatype<void(value_datatype_builder<T>&, T&)>::value
 || has_builder<T>::value
 >::type>
: std::true_type {};


// TODO: remove this after refactoring the building process for std::array, std::pair, std::tuple,
//       the custom struct macros and the builtin datatype


template <typename T>
inline typename std::enable_if<has_static_member_datatype<T, void(static_datatype_builder<T>&)>::value, datatype>::type
build_datatype() {
    //static_assert(!has_static_member_datatype<T, void(mxx::value_datatype_builder<T>&)>::value, "needs static datatype() function");
    //T val;
    mxx::static_datatype_builder<T> builder;
    T::datatype(builder);
    return builder.get_datatype();
}

template <typename T>
inline typename std::enable_if<
!has_static_member_datatype<T, void(static_datatype_builder<T>&)>::value
&& has_member_datatype<T, void(value_datatype_builder<T>&)>::value
, datatype>::type
build_datatype() {
    T val;
    value_datatype_builder<T> builder(val);
    val.datatype(builder);
    return builder.get_datatype();
}

// TODO: enable_if specializations for this function
template <typename T>
inline datatype build_datatype(const T&) {
    datatype dt(datatype_builder<T>::get_type(), is_builtin_type<T>::value);
    return dt;
}

// if datatype_builder<T> exists:
template <typename T>
inline typename std::enable_if<has_builder<T>::value, datatype>::type
build_datatype() {
    datatype dt(datatype_builder<T>::get_type(), is_builtin_type<T>::value);
    return dt;
}

template <typename T>
inline typename std::enable_if<!is_trivial_type<T>::value, datatype>::type
build_datatype() {
    // static assert the opposite to trigger the static assertion failure
    static_assert(is_trivial_type<T>::value,
    "Type `T` is not a `trivial` type and is thus not supported for mxx send/recv operations. "
    "This type needs one of the following to be supported as trivial datatype: "
    "specialized build_datatype<T>, a member function `datatype`, or global function `make_datatype(Layout& l, T&)`");

    return datatype();
}


template <typename T>
inline datatype get_datatype() {
    // TODO: retrieve cached datatype
    return build_datatype<T>();
}

template <typename T>
inline datatype get_datatype(const T& t) {
    // TODO: retrieve cached datatype
    return build_datatype(t);
}


/*********************************************************************
 *                      Custom struct datatypes                      *
 *********************************************************************/

// for each macros from: http://stackoverflow.com/questions/1872220/is-it-possible-to-iterate-over-arguments-in-variadic-macros
// Make a FOREACH macro
#define FE_1(WHAT, X) WHAT(X) 
#define FE_2(WHAT, X, ...) WHAT(X)FE_1(WHAT, __VA_ARGS__)
#define FE_3(WHAT, X, ...) WHAT(X)FE_2(WHAT, __VA_ARGS__)
#define FE_4(WHAT, X, ...) WHAT(X)FE_3(WHAT, __VA_ARGS__)
#define FE_5(WHAT, X, ...) WHAT(X)FE_4(WHAT, __VA_ARGS__)
#define FE_6(WHAT, X, ...) WHAT(X)FE_5(WHAT, __VA_ARGS__)
#define FE_7(WHAT, X, ...) WHAT(X)FE_6(WHAT, __VA_ARGS__)
#define FE_8(WHAT, X, ...) WHAT(X)FE_7(WHAT, __VA_ARGS__)
#define FE_9(WHAT, X, ...) WHAT(X)FE_8(WHAT, __VA_ARGS__)
#define FE_10(WHAT, X, ...) WHAT(X)FE_9(WHAT, __VA_ARGS__)

//... repeat as needed
#define GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,NAME,...) NAME 
#define FOR_EACH(action,...) \
  GET_MACRO(__VA_ARGS__,FE_10,FE_9,FE_8,FE_7,FE_6,FE_5,FE_4,FE_3,FE_2,FE_1)(action,__VA_ARGS__)

#define MXX_DT_PREAMBLE(BASE_TYPE) \
    MPI_Datatype _type; \
    BASE_TYPE p; \
    BASE_TYPE* pt = &p; \
    MPI_Aint p_adr; \
    MPI_Get_address(pt, &p_adr); \
    std::map<MPI_Aint, mxx::datatype> type_map;

#define MXX_DT_MEMBER_DISPLS(member) \
    MPI_Aint member ## _adr; \
    MPI_Get_address(&pt-> member, & member ## _adr); \
    type_map[member ## _adr - p_adr] = get_datatype<decltype(pt-> member )>();

#define MXX_DT_POSTAMBLE(BASE_TYPE) \
    int num_members = type_map.size(); \
    std::vector<int> blocklen(num_members, 1); \
    std::vector<MPI_Datatype> types(num_members); \
    std::vector<MPI_Aint> displs(num_members); \
    int i = 0; \
    for (auto& t : type_map) { \
      displs[i] = t.first; \
      types[i] = t.second.type(); \
      i++; \
    } \
    MPI_Datatype struct_type; \
    MPI_Type_create_struct((num_members), &blocklen[0], &displs[0], &types[0], &struct_type); \
    MPI_Type_create_resized(struct_type, 0, sizeof( BASE_TYPE ), &_type); \
    MPI_Type_commit(&_type); \
    MPI_Type_free(&struct_type); \
    return _type;


#define MXX_WRAP_TEMPLATE(...) __VA_ARGS__

#define MXX_DT_STRUCT_MEMBERS_GET_TYPE(BASE_TYPE, ...) MXX_DT_PREAMBLE(MXX_WRAP_TEMPLATE(BASE_TYPE)); FOR_EACH(MXX_DT_MEMBER_DISPLS, __VA_ARGS__); MXX_DT_POSTAMBLE(MXX_WRAP_TEMPLATE(BASE_TYPE));
#define MXX_DT_STRUCT_MEMBER_NUM_BASIC(MEMBER) datatype_builder<decltype(p. MEMBER)>::num_basic_elements()
#define MXX_DT_STRUCT_MEMBER_ADD_NUM_BASIC(MEMBER) + datatype_builder<decltype(p. MEMBER)>::num_basic_elements()

#define MXX_DT_STRUCT_MEMBERS_NUM_BASIC(BASE_TYPE, FIRST_MEMBER, ...) \
    static size_t num_basic_elements() { \
      BASE_TYPE p; \
      return MXX_DT_STRUCT_MEMBER_NUM_BASIC(FIRST_MEMBER) \
      FOR_EACH(MXX_DT_STRUCT_MEMBER_ADD_NUM_BASIC, __VA_ARGS__) ;\
    }



#define MXX_CUSTOM_STRUCT_(BASE_TYPE, ...) \
struct datatype_builder<BASE_TYPE> { \
    static MPI_Datatype get_type() { \
      MXX_DT_STRUCT_MEMBERS_GET_TYPE(MXX_WRAP_TEMPLATE(BASE_TYPE), __VA_ARGS__); \
    } \
    MXX_DT_STRUCT_MEMBERS_NUM_BASIC(MXX_WRAP_TEMPLATE(BASE_TYPE), __VA_ARGS__); \
};

#define MXX_CUSTOM_STRUCT(BASE_TYPE, ...) \
namespace mxx { \
template <> \
MXX_CUSTOM_STRUCT_(BASE_TYPE, __VA_ARGS__); \
} // namespace mxx


// use the MXX_WRAP_TEMPLATE() around templated types that have more than one paramter
// otherwise the comma "," in the template would split the templated type into separate arguments
#define MXX_CUSTOM_TEMPLATE_STRUCT(BASE_TYPE, ...) MXX_CUSTOM_STRUCT_(MXX_WRAP_TEMPLATE(BASE_TYPE), __VA_ARGS__)

} // namespace mxx


#endif // MXX_DATATYPES_HPP
