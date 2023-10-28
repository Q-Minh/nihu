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
 * \file relation.hpp 
 * \ingroup tmp
 * \brief Template metaprograms for comparison
 * \details
 * These comparisons are used in comparing asymptotic singularity types, e.g.
 */ 

#ifndef RELATION_HPP_INCLUDED
#define RELATION_HPP_INCLUDED

#include <type_traits>

namespace tmp
{
	/** 
	 * \brief General declaration of the less oparation 
	 * \tparam N Left hand side of comparison
	 * \tparam M Right hand side of comparison
	 * \details 
	 * Return \c std::true_type if first is less than second, i.e., \c N < \c M
	 */
	template <class N, class M> struct less;

	/** 
	 * \brief General declaration of the greater oparation 
	 * \tparam N Left hand side of comparison
	 * \tparam M Right hand side of comparison
	 * \details 
	 * Return \c std::true_type if first is greater than second, i.e.,
	 * \c N > \c M
	 */
	template <class N, class M> struct greater;

	/** 
	 * \brief Compute maximum of types
	 * \tparam Val First value in the list
	 * \tparam Args Other values in the list
	 * \details
	 * This is the general recursive branch of the template metaprogram.
	 * The recursion is implemented by taking the maximum of first value and the
	 * maximum of the other values in the list.
	 */
	template <class Val, class...Args>
	struct max_ : max_<Val, typename max_<Args...>::type> {};

	/** 
	 * \brief Compute maximum of type, specialisation for two parameters
	 * \tparam Val1 Left hand side of the comparison
	 * \tparam Val2 Right hand side of the comparison
	 * \details
	 * Returns the greater (see \ref greater) of the two parameters.
	 */
	template <class Val1, class Val2>
	struct max_<Val1, Val2> : std::conditional<
		greater<Val1, Val2>::value,
		Val1, Val2
	> {};

	/** 
	 * \brief Compute minimum of types 
	 * \tparam Val First value in the list
	 * \tparam Args Other values in the list
	 * \details
	 * This is the general recursive branch of the template metaprogram.
	 * The recursion is implemented by taking the minimum of first value and the
	 * minimum of the other values in the list.
	 */
	template <class Val, class...Args>
	struct min_ : min_<Val, typename min_<Args...>::type> {};

	/** 
	 * \brief Compute minimum of types, specialisation for two parameters 
	 * \tparam Val1 Left hand side of comparison
	 * \tparam Val2 Right hand side of comparison
	 * \details
	 * Returns the smaller (see \ref less) of the two parameters.
	 */
	template <class Val1, class Val2>
	struct min_<Val1, Val2> : std::conditional<
		less<Val1, Val2>::value,
		Val1, Val2
	> {};
}

#endif /* RELATION_HPP_INCLUDED */
