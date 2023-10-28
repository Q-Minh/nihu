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
 * \brief Plain type calculations
 * \ingroup util
 */

#ifndef PLAIN_TYPE_HPP_INCLUDED
#define PLAIN_TYPE_HPP_INCLUDED

#include "eigen_utils.hpp"

#include <type_traits>

namespace NiHu
{

/** 
 * \brief Plain object type of a class
 * \tparam T the class to convert to plain type
 * \details
 * This is the general case where the class is not an Eigen expression
 */
template <
	class T,
	bool isEigen = is_eigen<T>::value
>
struct plain_type 
	: std::decay<T>
{
};


/** 
 * \brief Plain object type of a class
 * \tparam T the expression class to convert to plain type
 * \details
 * Specialisation for the case of Eigen expressions
 */
template <class T>
struct plain_type<T, true>
{
	typedef typename std::decay<T>::type::PlainObject type;
};

} // end of namespace NiHu

#endif /* PLAIN_TYPE_HPP_INCLUDED */

