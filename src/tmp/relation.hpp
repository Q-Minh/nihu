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

#ifndef RELATION_HPP_INCLUDED
#define RELATION_HPP_INCLUDED

#include <type_traits>

namespace tmp
{
	/** \brief return true_type if first is less than second */
	template <class N, class M> struct less;

	/** \brief return true_type if first is greater than second */
	template <class N, class M> struct greater;


	/** \brief compute maximum of types */
	template <class Val, class...Args>
	struct max_ : max_<Val, typename max_<Args...>::type> {};

	/** \brief specialisation of max_ for the two parameter case */
	template <class Val1, class Val2>
	struct max_<Val1, Val2> : std::conditional<
		greater<Val1, Val2>::value,
		Val1, Val2
	> {};

	/** \brief compute minimum of types */
	template <class Val, class...Args>
	struct min_ : min_<Val, typename min_<Args...>::type> {};

	/** \brief specialisation of min_ for the two parameter case */
	template <class Val1, class Val2>
	struct min_<Val1, Val2> : std::conditional<
		less<Val1, Val2>::value,
		Val1, Val2
	> {};
}


#endif // RELATION_HPP_INCLUDED
