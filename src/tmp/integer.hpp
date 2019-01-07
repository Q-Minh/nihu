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
 * \file integer.hpp
 * \brief integer type representation and basic integer arithmetics
 */

#ifndef INTEGER_HPP
#define INTEGER_HPP

#include "relation.hpp"
#include "operator.hpp"
//#include "bool.hpp"

namespace tmp
{
	/**
	 * \brief integer type representation
	 */
	template <class T, T N>
	struct integer : std::integral_constant<T, N>
	{
		/** \brief self returning metafunction */
		typedef integer type;
		/** \brief the next value */
		static T const next = N+1;
		/** \brief the previous value */
		static T const prev = N-1;
	};

	/**
	 * \brief metafunction returning next integer
	 */
	template <class T, T N>
	struct next<integer<T, N> > : integer<T, integer<T, N>::next> {};

	/**
	 * \brief metafunction returning previous integer
	 */
	template <class T, T N>
	struct prev<integer<T, N> > : integer<T, integer<T, N>::prev> {};

	/**
	 * \brief metafunction returning sum of two integers
	 */
	template <class T, T N, T M>
	struct plus<integer<T, N>, integer<T, M> > : integer<T, N+M> {};

	/**
	 * \brief metafunction returning difference of two integers
	 */
	template <class T, T N, T M>
	struct minus<integer<T, N>, integer<T, M> > : integer<T, N-M> {};

	/**
	 * \brief metafunction returning difference of two integers
	 */
	template <class T, T N, T M>
	struct mul<integer<T, N>, integer<T, M> > : integer<T, N*M> {};

	template <class T, T N, T M>
	struct less<integer<T, N>, integer<T, M> > : std::integral_constant<bool, (N<M) > {};

	template <class T, T N, T M>
	struct greater<integer<T, N>, integer<T, M> > : std::integral_constant<bool, (N>M) > {};
}

#endif

