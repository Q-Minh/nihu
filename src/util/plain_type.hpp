/**
 * \file plain_type.hpp
 * \brief plain type calculations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef PLAIN_TYPE_HPP_INCLUDED
#define PLAIN_TYPE_HPP_INCLUDED

#include <type_traits>
#include "../bem/couple.hpp"
#include "Eigen/Dense"

template <class T>
struct is_eigen : std::is_base_of<Eigen::EigenBase<T>, T> {};

template <class T, bool isEigen = is_eigen<typename std::decay<T>::type>::value, bool isCouple = std::is_base_of<couple_base<T>, T>::value>
struct plain_type
{
	typedef T type;
};

template <class T>
struct plain_type<T, true, false> : std::decay<
	decltype(static_cast<
		typename std::decay<T>::type *
	>(nullptr)->eval())
>{};

template <class T>
struct plain_type<T, false, true>
{
	typedef couple<
		typename plain_type<
			typename T::first_t
		>::type,
		typename plain_type<
			typename T::second_t
		>::type
	> type;
};

#endif // PLAIN_TYPE_HPP_INCLUDED

