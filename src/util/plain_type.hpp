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
template <class T, bool isEigen = is_eigen<T>::value, bool isCouple = is_couple<T>::value>
struct plain_type : std::decay<T> {};


/** \brief specialisation of ::plain_type for the case of eigen expressions
 * \tparam T the expression class to convert to plain type
 */
template <class T>
struct plain_type<T, true, false>
{
	typedef typename T::PlainObject type;
};

/** \brief specialisation of ::plain_type for the case of a couple expression
 * \tparam T the couple expression class to convert to plain type
 * \todo should work recursively and for couple expressions, not only for couples
 */
template <class T>
struct plain_type<T, false, true> : couple<
	typename plain_type<
		typename T::template couple_type<0>::type
	>::type,
	typename plain_type<
		typename T::template couple_type<1>::type
	>::type
> {};

#endif // PLAIN_TYPE_HPP_INCLUDED

