#ifndef MXX_TYPE_TRAITS
#define MXX_TYPE_TRAITS

#include <type_traits>
#include <vector>
#include <string>

namespace mxx {

// source for testing if member functions are available: http://stackoverflow.com/a/16824239/4639394
#define MXX_DEFINE_HAS_MEMBER(member) \
template<typename, typename T> \
struct has_member_ ## member {}; \
\
template<typename C, typename Ret, typename... Args> \
struct has_member_ ## member <C, Ret(Args...)> { \
private: \
    template<typename T> \
    static constexpr auto check(T*) \
    -> typename \
        std::is_same< \
            decltype( std::declval<T>(). member ( std::declval<Args>()... ) ), \
            Ret    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */  \
        >::type;  /* attempt to call it and see if the return type is correct */  \
    template <typename> \
    static constexpr std::false_type check(...); \
    typedef decltype(check<C>(0)) type; \
public: \
    static constexpr bool value = type::value; \
};
#define MXX_DEFINE_HAS_STATIC_MEMBER(member) \
template<typename, typename T> \
struct has_static_member_ ## member {}; \
\
template<typename C, typename Ret, typename... Args> \
struct has_static_member_ ## member <C, Ret(Args...)> { \
private: \
    template<typename T> \
    static constexpr auto check(T*) \
    -> typename \
        std::is_same< \
            decltype(T:: member ( std::declval<Args>()... ) ), \
            Ret    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */  \
        >::type;  /* attempt to call it and see if the return type is correct */  \
    template<typename> \
    static constexpr std::false_type check(...); \
    typedef decltype(check<C>(0)) type; \
public: \
    static constexpr bool value = type::value; \
};

#define MXX_DEFINE_IS_GLOBAL_FUNC(fname) \
template<typename U> \
struct is_global_func_ ## fname { \
private: \
    typedef U signature; \
    template <typename T, T> struct has_matching_sig; \
    template <typename T> \
    static std::true_type check(has_matching_sig<T*,& :: fname>*); \
    template <typename T> \
    static std::false_type check(...); \
    typedef decltype(check<signature>(0)) type; \
public: \
    static constexpr bool value = type::value; \
};

// own implementation of the HAS_MEMBER struct
template <typename, typename>
struct has_member_size : std::false_type {};

template <typename C, typename R, typename... Args>
struct has_member_size<C, R(Args...)> : std::integral_constant<bool,
    std::is_same<decltype(std::declval<C>().size(std::declval<Args>()...)), R>::value> {};

// type traits returning whether the type `C` has the typedef member ::value_type
#define MXX_DEFINE_HAS_TYPEDEF(type_name) \
template <typename C, typename Enable = void> \
struct has_typedef_ ## type_name : std::false_type {}; \
template <typename C> \
struct has_typedef_ ## type_name <C, typename std::enable_if< \
!std::is_same<typename C:: type_name ,void>::value>::type> \
: std::true_type {};

MXX_DEFINE_HAS_TYPEDEF(value_type)
MXX_DEFINE_HAS_TYPEDEF(iterator)
MXX_DEFINE_HAS_TYPEDEF(const_iterator)

MXX_DEFINE_HAS_MEMBER(data)
//MXX_DEFINE_HAS_MEMBER(size)
MXX_DEFINE_HAS_MEMBER(resize)

MXX_DEFINE_HAS_MEMBER(end)
MXX_DEFINE_HAS_MEMBER(begin)

MXX_DEFINE_HAS_STATIC_MEMBER(datatype)
MXX_DEFINE_HAS_MEMBER(datatype)


// TODO: build this into the mxx/datatype structures
template <typename T, typename Enable = void>
struct has_datatype : std::false_type {};



/***********************************************************
 *  Compile-time reflection for members begin() and end()  *
 ***********************************************************/

template <typename C, typename Enable = void>
struct has_nonconst_begin_end : std::false_type {};

template <typename C>
struct has_nonconst_begin_end<C, typename std::enable_if<has_typedef_iterator<C>::value>::type>
: std::integral_constant<bool,
    has_member_begin<typename std::remove_const<C>::type, typename C::iterator()>::value
    && has_member_end<typename std::remove_const<C>::type, typename C::iterator()>::value> {};

template <typename C, typename Enable = void>
struct has_const_begin_end : std::false_type {};

template <typename C>
struct has_const_begin_end<C, typename std::enable_if<has_typedef_iterator<C>::value>::type>
: std::integral_constant<bool,
    has_member_begin<const typename std::remove_const<C>::type, typename C::const_iterator()>::value
    && has_member_end<const typename std::remove_const<C>::type, typename C::const_iterator()>::value> {};

template <typename C>
struct has_begin_end : std::integral_constant<bool,
    has_nonconst_begin_end<C>::value
    && has_const_begin_end<C>::value> {};

/**************************************************************
 *  Compile-time determination whether a type is a container  *
 **************************************************************/

/*
template <typename C, typename Enable = void>
struct is_container : std::false_type {};
*/
template <typename C, typename Enable = void>
struct is_container : std::integral_constant<bool,
    has_begin_end<C>::value> {};

template <typename C, typename Enable = void>
struct is_flat_container : std::false_type {};

template <typename C>
struct is_flat_container<C, typename std::enable_if<has_typedef_value_type<C>::value>::type>
: std::integral_constant<bool,
    is_container<C>::value
    && mxx::has_datatype<typename C::value_type>::value> {};

template <typename C, typename Enable = void>
struct is_contiguous_container : std::false_type {};

template <typename T, typename A>
struct is_contiguous_container<std::vector<T,A>> : std::true_type {};

template <typename C, typename CT, typename A>
struct is_contiguous_container<std::basic_string<C, CT, A>> : std::true_type {};


// anything that's templated by at least its value type T and has the member functions
// size_t size(), void resize(size_t size), and T* data()
template <template <typename, typename...> class C, typename T, typename... Targs>
struct is_contiguous_container<C<T,Targs...>> : std::integral_constant<bool,
        has_member_data<C<T,Targs...>, T*()>::value &&
        has_member_size<C<T,Targs...>, std::size_t()>::value &&
        has_member_resize<C<T,Targs...>, void(std::size_t)>::value> {};


template <typename C>
struct is_flat_contiguous_container : std::integral_constant<bool,
    is_contiguous_container<C>::value && is_flat_container<C>::value> {};

template <typename C, typename Enable = void>
struct is_flat_type : std::integral_constant<bool,
    is_contiguous_container<C>::value> {};

} // namespace mxx

#endif // MXX_TYPE_TRAITS
