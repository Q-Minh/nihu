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
 * \file lib_element.hpp 
 * \ingroup lib_sef
 */

#ifndef LIB_ELEMENT_HPP_INCLUDED
#define LIB_ELEMENT_HPP_INCLUDED

#include "../core/element.hpp"
#include "../util/type2tag.hpp"
#include "lib_shape.hpp"

namespace NiHu
{

/** \brief A linear (2-noded) line element in 2D space */
typedef surface_element<line_1_shape_set, space_2d<>::scalar_t> line_1_elem;

/** \brief A quadratic (3-noded) line element in 2D space */
typedef surface_element<line_2_shape_set, space_2d<>::scalar_t> line_2_elem;

/** \brief A linear (3-noded) triangular element in 3D space */
typedef surface_element<tria_1_shape_set, space_3d<>::scalar_t> tria_1_elem;

/** \brief A quadratic (6-noded) triangular element in 3D space */
typedef surface_element<tria_2_shape_set, space_3d<>::scalar_t> tria_2_elem;

/** \brief A linear (4-noded) quadrangular element in 3D space */
typedef surface_element<quad_1_shape_set, space_3d<>::scalar_t> quad_1_elem;

/** \brief A quadratic (9-noded) quadrangular element in 3D space */
typedef surface_element<quad_2_shape_set, space_3d<>::scalar_t> quad_2_elem;

/** \brief A quadratic (8-noded) quadrangular element in 3D space */
typedef surface_element<quad_28_shape_set, space_3d<>::scalar_t> quad_28_elem;


/** \brief A linear (2-noded) linear volume element in 1D space */
typedef volume_element<line_1_shape_set, space_1d<>::scalar_t> line_1_volume_elem;

/** \brief A linear (3-noded) triangular volume element in 2D space */
typedef volume_element<tria_1_shape_set, space_2d<>::scalar_t> tria_1_volume_elem;

/** \brief A linear (4-noded) quadrilateral volume element in 2D space */
typedef volume_element<quad_1_shape_set, space_2d<>::scalar_t> quad_1_volume_elem;

typedef type2tag<line_1_elem> line_1_tag;
typedef type2tag<line_2_elem> line_2_tag;
typedef type2tag<tria_1_elem> tria_1_tag;
typedef type2tag<tria_2_elem> tria_2_tag;
typedef type2tag<quad_1_elem> quad_1_tag;
typedef type2tag<quad_2_elem> quad_2_tag;
typedef type2tag<quad_28_elem> quad_28_tag;
typedef type2tag<line_1_volume_elem> line_1_volume_tag;
typedef type2tag<tria_1_volume_elem> tria_1_volume_tag;
typedef type2tag<quad_1_volume_elem> quad_1_volume_tag;

}

#endif // LIB_ELEMENT_HPP_INCLUDED

