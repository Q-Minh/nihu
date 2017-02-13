// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file plain_type.hpp
 * \brief plain type calculations
 */

#ifndef PLAIN_TYPE_HPP_INCLUDED
#define PLAIN_TYPE_HPP_INCLUDED

#include <type_traits>
#include "couple.hpp"
#include "eigen_utils.hpp"

namespace NiHu
{

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

}

#endif // PLAIN_TYPE_HPP_INCLUDED

