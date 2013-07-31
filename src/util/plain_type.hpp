/**
 * \file plain_type.hpp
 * \brief plain type calculations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef PLAIN_TYPE_HPP_INCLUDED
#define PLAIN_TYPE_HPP_INCLUDED

#include <type_traits>
#include "couple.hpp"
#include "eigen_utils.hpp"

/** \brief plain object type of a class
 * \tparam T the class to convert to plain type
 * \details this is the general case where the class is not an expression
 */
template <
	class T,
	bool isEigen = is_eigen<typename std::decay<T>::type>::value,
	bool isCouple = is_couple<typename std::decay<T>::type>::value
>
struct plain_type : std::decay<T> {};


/** \brief specialisation of ::plain_type for the case of eigen expressions
 * \tparam T the expression class to convert to plain type
 */
template <class T>
struct plain_type<T, true, false>
{
	typedef typename std::decay<T>::type::PlainObject type;
};


template <class T>
struct tuple_plain;

template <class...Args>
struct tuple_plain<std::tuple<Args...> > : couple<
	typename plain_type<Args>::type...
> {};

/** \brief specialisation of ::plain_type for the case of a couple expression
 * \tparam T the couple expression class to convert to plain type
 */
template <class T>
struct plain_type<T, false, true> : tuple_plain<
	typename couple_traits<T>::tuple_t
> {};
/*
template <class T>
struct plain_type<T, false, true> : couple<
	typename plain_type<
		decltype( static_cast<typename std::decay<T>::type const *>(nullptr)->template get<0>() )
	>::type,
	typename plain_type<
		decltype( static_cast<typename std::decay<T>::type const *>(nullptr)->template get<1>() )
	>::type
> {};
*/

#endif // PLAIN_TYPE_HPP_INCLUDED

