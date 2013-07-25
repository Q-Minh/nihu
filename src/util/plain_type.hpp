/**
 * \file plain_type.hpp
 * \brief plain type calculations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef PLAIN_TYPE_HPP_INCLUDED
#define PLAIN_TYPE_HPP_INCLUDED

#include <type_traits>
#include "couple.hpp"
#include "Eigen/Dense"

/** \brief metafunction determining if its argument is an Eigen expression or not
 * \tparam T the class to investigate
 */
template <class T>
struct is_eigen : std::is_base_of<
	Eigen::EigenBase<typename std::decay<T>::type>,
	typename std::decay<T>::type
> {};

template <class m1, class m2, int t>
struct is_eigen<Eigen::GeneralProduct<m1, m2, t> > : std::true_type {};

/** \brief plain object type of a class
 * \tparam T the class to convert to plain type
 * \details this is the general case where the class is not an expression
 */
template <class T, bool isEigen = is_eigen<T>::value, bool isCouple = is_couple<T>::value>
struct plain_type : std::decay<T> {};

/** \brief specialisation of ::plain_type for the case of an eigen expression
 * \tparam T the expression class to convert to plain type
 * \details plain_type is the result type of function eval()
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

