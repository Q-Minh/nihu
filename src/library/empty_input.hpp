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
 * \file empty_input.hpp
 * \ingroup library
 * \brief implementation of empty kernel
 */

#ifndef EMPTY_INPUT_HPP_INCLUDED
#define EMPTY_INPUT_HPP_INCLUDED

namespace NiHu
{

// forward declaration
template <class space_t>
class empty_input;

/** \brief traits of a location
* \tparam Space the x-space
*/
template <class Space>
struct kernel_input_traits<empty_input<Space> >
{
	/** \brief the space type */
	typedef Space space_t;
};


/** \brief a kernel input representing a single location \f$ \gamma = {\bf x} \f$
* \tparam Space the coordinate space of the location
*/
template <class Space>
class empty_input :
	public kernel_input_base<empty_input<Space> >
{
public:
	/** \brief the CRTP base class */
	typedef kernel_input_base<empty_input<Space> > base_t;

	/** \brief constructor from element and reference domain vector
	* \tparam elem_t the element type
	* \param [in] elem the element to construct from
	* \param [in] xi the location in the reference domain
	*/
	template <class elem_t>
	empty_input(elem_t const &elem, typename elem_t::xi_t const &xi)
		: base_t(elem, xi)
	{
	}
};

}


#endif // EMPTY_INPUT_HPP_INCLUDED

