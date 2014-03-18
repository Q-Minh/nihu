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
 * \file operator.hpp
 * \brief declaration of operator metafunctions
 * \details This file only declares the binary operators so that they can be overloaded in specific applications
 */
 
#ifndef OPERATOR_HPP
#define OPERATOR_HPP

namespace tmp
{
	/** \brief increment operator */
	template <class A>
	struct next;

	/** \brief decrement operator */
	template <class A>
	struct prev;

	/** \brief binary plus */
	template <class A, class B>
	struct plus;

	/** \brief binary minus */
	template <class A, class B>
	struct minus;

	/** \brief binary multiply */
	template <class A, class B>
	struct mul;

	/** \brief binary less than */
	template <class A, class B>
	struct less;
}

#endif

