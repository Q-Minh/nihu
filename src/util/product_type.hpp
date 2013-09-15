// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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
 * \file product_type.hpp
 * \brief product type calculations
 */

#ifndef PRODUCT_TYPE_HPP_INCLUDED
#define PRODUCT_TYPE_HPP_INCLUDED

#include "eigen_utils.hpp"

/** \brief metafunction returning the product type of two classes
 * \tparam Lhs the left hand side expression type
 * \tparam Rhs the right hand side expression type
 */
template<class Lhs, class Rhs, bool isBothEigen = is_eigen<Lhs>::value && is_eigen<Rhs>::value>
struct product_type
{
	/** \brief the return type computed by decltype */
	typedef decltype(
		(*static_cast<typename std::decay<Lhs>::type *>(nullptr))
		*
		(*static_cast<typename std::decay<Rhs>::type *>(nullptr))
	) type;
};


/** \brief specialisation of ::product_type for the case of two eigen expressions
 * \tparam Lhs the left hand side eigen expression type
 * \tparam Rhs the right hand side eigen expression type
 */
template<class Lhs, class Rhs>
struct product_type<Lhs, Rhs, true>
{
	typedef typename Eigen::ProductReturnType<Lhs, Rhs>::Type type;
};

#endif // PRODUCT_TYPE_HPP_INCLUDED

