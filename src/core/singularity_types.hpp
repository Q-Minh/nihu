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

/** \file singularity_types.hpp
 * \brief definition of different singularity types
 */

#ifndef SINGULARITY_TYPES_HPP_INCLUDED
#define SINGULARITY_TYPES_HPP_INCLUDED

#include "../tmp/relation.hpp"

/** \brief namespace encapsulating singularity type classes */
namespace singularity_type
{
	/** \brief no singularity */
	struct regular {};

	/** \brief logarithmic singularity with specified order \f$ log(r)^o \f$
	 * \tparam order the order of the logarithm
	 */
	template <unsigned order>
	struct log {};

	/** \brief inverse singularity with specified order \f$ 1/r^o \f$
	 * \tparam order the singularity order
	 */
	template <unsigned order>
	struct inverse {};
}

/** \brief returns the minimal reference domain dimension where the singularity can be integrated
 * \tparam SingularityType the singularity type
 */
template <class SingularityType>
struct minimal_reference_dimension;

/** \brief specialisation of ::minimal_reference_dimension to the log<1> singularity */
template <>
struct minimal_reference_dimension<singularity_type::log<1> >
{
	static unsigned const value = 1;
};


/** \brief specialisation of ::minimal_reference_dimension to the log<o> singularity
 * \tparam order the inverse singularity order
 */
template <unsigned order>
struct minimal_reference_dimension<singularity_type::inverse<order> >
{
	static unsigned const value = order+1;
};



namespace tmp
{
	template <unsigned o1, unsigned o2>
	struct less<singularity_type::log<o1>, singularity_type::log<o2> >
		: std::integral_constant<bool, (o1 < o2) > {};

	template <unsigned o1, unsigned o2>
	struct greater<singularity_type::log<o1>, singularity_type::log<o2> >
		: std::integral_constant<bool, (o1 > o2) > {};

	template <unsigned o1, unsigned o2>
	struct less<singularity_type::inverse<o1>, singularity_type::inverse<o2> >
		: std::integral_constant<bool, (o1 < o2) > {};

	template <unsigned o1, unsigned o2>
	struct greater<singularity_type::inverse<o1>, singularity_type::inverse<o2> >
		: std::integral_constant<bool, (o1 > o2) > {};

	template <unsigned o1, unsigned o2>
	struct less<singularity_type::log<o1>, singularity_type::inverse<o2> >
		: std::true_type {};

	template <unsigned o1, unsigned o2>
	struct less<singularity_type::inverse<o1>, singularity_type::log<o2> >
		: std::false_type {};

	template <unsigned o1, unsigned o2>
	struct greater<singularity_type::log<o1>, singularity_type::inverse<o2> >
		: std::false_type {};

	template <unsigned o1, unsigned o2>
	struct greater<singularity_type::inverse<o1>, singularity_type::log<o2> >
		: std::true_type {};
}


#endif // SINGULARITY_TYPES_HPP_INCLUDED
