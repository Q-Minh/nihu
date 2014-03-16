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

#include "lib_shape.hpp"

namespace shape_set_traits
{
	template <>
	const std::string name<line_0_shape_set>::value = "Line 0 shape set";

	template <>
	const std::string name<line_1_shape_set>::value = "Line 1 shape set";

	template <>
	const std::string name<line_2_shape_set>::value = "Line 2 shape set";

	template <>
	const std::string name<tria_0_shape_set>::value = "Tria 0 shape set";

	template <>
	const std::string name<tria_1_shape_set>::value = "Tria 1 shape set";

	template <>
	const std::string name<tria_2_shape_set>::value = "Tria 2 shape set";

	template <>
	const std::string name<quad_0_shape_set>::value = "Quad 0 shape set";

	template <>
	const std::string name<quad_1_shape_set>::value = "Quad 1 shape set";

	template <>
	const std::string name<quad_2_shape_set>::value = "Quad 2 shape set";

	template <>
	const std::string name<quad_28_shape_set>::value = "Quad 28 shape set";

	template <>
	const std::string name<brick_0_shape_set>::value = "Brick 0 shape set";

	template <>
	const std::string name<brick_1_shape_set>::value = "Brick 1 shape set";

	template <>
	const std::string name<parallelogram_shape_set>::value = "Parallelogram shape set";
}


line_2_shape_set::xi_t
const line_2_shape_set::m_corners[line_2_shape_set::num_nodes] = {
	line_2_shape_set::xi_t::Constant(-1.0),
	line_2_shape_set::xi_t::Constant(0.0),
	line_2_shape_set::xi_t::Constant(1.0),
};

int const line_2_shape_set::m_domain_indices[line_2_shape_set::num_nodes] = { 0, -1, 1 };


tria_2_shape_set::xi_t
const tria_2_shape_set::m_corners[tria_2_shape_set::num_nodes] = {
	tria_2_shape_set::xi_t(0.0, 0.0),
	tria_2_shape_set::xi_t(0.5, 0.0),
	tria_2_shape_set::xi_t(1.0, 0.0),
	tria_2_shape_set::xi_t(0.5, 0.5),
	tria_2_shape_set::xi_t(0.0, 1.0),
	tria_2_shape_set::xi_t(0.0, 0.5)
};

int const tria_2_shape_set::m_domain_corners[tria_2_shape_set::num_nodes] =
{ 0, -1, 1, -1, 2, -1 };


quad_2_shape_set::xi_t
const quad_2_shape_set::m_corners[quad_2_shape_set::num_nodes] = {
	quad_2_shape_set::xi_t(-1.0, -1.0),
	quad_2_shape_set::xi_t(0.0, -1.0),
	quad_2_shape_set::xi_t(+1.0, -1.0),
	quad_2_shape_set::xi_t(+1.0, 0.0),
	quad_2_shape_set::xi_t(+1.0, +1.0),
	quad_2_shape_set::xi_t(0.0, +1.0),
	quad_2_shape_set::xi_t(-1.0, +1.0),
	quad_2_shape_set::xi_t(-1.0, 0.0),
	quad_2_shape_set::xi_t(0.0, 0.0)
};

int const quad_2_shape_set::m_domain_corners[quad_2_shape_set::num_nodes] =
{ 0, -1, 1, -1, 2, -1, 3, -1, -1 };


quad_28_shape_set::xi_t
const quad_28_shape_set::m_corners[quad_28_shape_set::num_nodes] = {
	quad_28_shape_set::xi_t(-1.0, -1.0),
	quad_28_shape_set::xi_t(0.0, -1.0),
	quad_28_shape_set::xi_t(+1.0, -1.0),
	quad_28_shape_set::xi_t(+1.0, 0.0),
	quad_28_shape_set::xi_t(+1.0, +1.0),
	quad_28_shape_set::xi_t(0.0, +1.0),
	quad_28_shape_set::xi_t(-1.0, +1.0),
	quad_28_shape_set::xi_t(-1.0, 0.0)
};

int const quad_28_shape_set::m_domain_corners[quad_28_shape_set::num_nodes] =
{ 0, -1, 1, -1, 2, -1, 3, -1 };

