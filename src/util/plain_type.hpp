/**
 * \file plain_type.hpp
 * \brief plain type calculations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef PLAIN_TYPE_HPP_INCLUDED
#define PLAIN_TYPE_HPP_INCLUDED

#include <type_traits>
#include "../bem/couple.hpp"
#include "product_type.hpp"

template <class T, bool isEigen = is_eigen_matrix<T>::value, bool isCouple = is_couple<T>::value>
struct plain_type
{
	typedef T type;
};


template <class T>
struct plain_type<T, true, false>
{
	typedef typename T::PlainObject type;
};


template <class T>
struct plain_type<T, false, true>
{
	typedef couple<
		typename plain_type<typename T::first_t>::type,
		typename plain_type<typename T::second_t>::type
	> type;
};


#endif // PLAIN_TYPE_HPP_INCLUDED
