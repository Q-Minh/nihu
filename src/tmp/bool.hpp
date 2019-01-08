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
 * \file bool.hpp
 * \ingroup tmp
 * \brief implementation of Boolean functions
 * \details 
 * The boolean type is already contained in the standard \c type_traits library 
 * as std::integral_constant<bool, x>.
 * 
 * This file implements some Boolean functions and the compile time \c if_ 
 * control structure.
 */

#ifndef BOOL_HPP_INCLUDED
#define BOOL_HPP_INCLUDED

#include <type_traits>

namespace tmp
{
	/**
	 * \brief Boolean negation
	 * \tparam A Boolean data
	 * \returns not A
	 */
	template <class A>
	struct not_ : std::integral_constant<bool, !A::value> {};

	/** 
	 * \brief Disjunction of Boolean constants 
	 * \tparam Args Boolean values to disjunct
	 * \returns Arg1 OR Arg2 OR ... OR ArgN 
	 * \details
	 * This is the general case
	 */
	template <class...Args>
	struct or_ : std::false_type {};
	
	/**
	 * \brief Disjunction of Boolean constants
	 * \tparam Args Booen values to disjunct
	 * \returns Arg1 OR Arg2 OR ... OR ArgN
	 * \details 
	 * This partial specialisation covers the case False OR Arg2 ...
	 */
	template <class...Args>
	struct or_<std::false_type, Args...> : or_<Args...> {};

	/**
	 * \brief Disjunction of Boolean constants
	 * \tparam Args Booen values to disjunct
	 * \returns Arg1 OR Arg2 OR ... OR ArgN
	 * \details 
	 * This partial specialisation covers the case True OR Arg2 ...
	 */
	template <class...Args>
	struct or_<std::true_type, Args...> : std::true_type {};

	/** 
	 * \brief Conjunction of boolean constants
	 * \tparam Args Boolen values to conjugate
	 * \returns Arg1 AND Arg2 AND ... AND ArgN
	 * \details 
	 * This is the general case
	 */ 
	template <class...Args>
	struct and_ : std::false_type {};

	/** 
	 * \brief Conjunction of boolean constants (recursion terminator)
	 * \returns True
	 * \details 
	 * This partial specialisation covers the case True to terminate the 
	 * the recursion used in the metafunction. 
	 */ 
	template <>
	struct and_<std::true_type> : std::true_type {};
	
	/** 
	 * \brief Conjunction of boolean constants
	 * \tparam Args Boolen values to conjugate
	 * \returns Arg1 AND Arg2 AND ... AND ArgN
	 * \details 
	 * This partial specialisation covers the case True AND Arg2 ...
	 */ 
	template <class...Args>
	struct and_<std::true_type, Args...> : and_<Args...> {};

	/** 
	 * \brief Conjunction of boolean constants
	 * \tparam Args Boolen values to conjugate
	 * \returns Arg1 AND Arg2 AND ... AND ArgN
	 * \details 
	 * This partial specialisation covers the case False AND Arg2 ...
	 */ 
	template <class...Args>
	struct and_<std::false_type, Args...> : std::false_type {};
	
}

#endif /* BOOL_HPP_INCLUDED */

