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

#ifndef QUAD_28_ELEM_HPP_INCLUDED
#define QUAD_28_ELEM_HPP_INCLUDED

#include "quad_28_shape_set.hpp"
#include "../core/element.hpp"

/** \brief quadratic 9-noded quadrilateral element */
typedef general_surface_element<quad_28_shape_set, quad_1_elem::space_t::scalar_t> quad_28_elem;

/** \brief tag of a 9-noded quadratic quad element */
struct quad_28_tag {};

/** \brief element assigned to the quad_28_tag */
template <>
struct tag2element<quad_28_tag>
{
	typedef quad_28_elem type;
};

#endif
