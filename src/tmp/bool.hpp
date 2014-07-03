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
 * \details The boolean type is already contained in type_traits as std::integral_constant<bool, x>.
 * This file implements some Boolean functions and the compile time if_ control structure
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

	/** \brief disjunction of boolean constants */
	template <class...Args>
	struct or_ : std::false_type {};

	template <class...Args>
	struct or_<std::false_type, Args...> : or_<Args...> {};

	template <class...Args>
	struct or_<std::true_type, Args...> : std::true_type {};

	/** \brief conjunction of boolean constants */
	template <class...Args>
	struct and_ : std::false_type {};

	template <class...Args>
	struct and_<std::true_type, Args...> : and_<Args...> {};

	template <class...Args>
	struct and_<std::false_type, Args...> : std::false_type {};

	template <>
	struct and_<std::true_type> : std::true_type {};

	/**
	 * \brief IF control structure
	 * \tparam Choice a choice evaluated to a logical type
	 * \tparam T the type returned when Choice is true_type
	 * \tparam F the type returned when Choice is false_type
	 */
	template <class Cond, class T, class F>
	struct if_ : std::conditional<Cond::value, T, F> {};
}

#endif

