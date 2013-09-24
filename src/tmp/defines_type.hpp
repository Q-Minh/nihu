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
 * \file defines_type.hpp
 * \brief metafunction determining typedefs
 */

#ifndef DEFINES_TYPE_HPP_INCLUDED
#define DEFINES_TYPE_HPP_INCLUDED

#include <type_traits>

namespace internal
{
	/** \brief metafunction converting anything to void */
	template <class C>
	struct voidize { typedef void type; };
}

#define DEFINES_TYPE(TYPE)								\
template <class C, class = void>						\
struct defines_##TYPE : std::false_type {};				\
template <class C>										\
														\
struct defines_##TYPE <									\
	C,													\
	typename internal::voidize<typename C::TYPE>::type	\
> : std::true_type { };

#endif // DEFINES_TYPE_HPP_INCLUDED
