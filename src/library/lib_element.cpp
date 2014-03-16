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

#include "lib_element.hpp"

namespace element_traits
{
	template <>
	const std::string name<line_1_elem>::value = "Line 1 elem";
	template <>
	const std::string name<line_2_elem>::value = "Line 2 elem";
	template <>
	const std::string name<tria_1_elem>::value = "Tria 1 elem";
	template <>
	const std::string name<tria_2_elem>::value = "Tria 2 elem";
	template <>
	const std::string name<quad_1_elem>::value = "Quad 1 elem";
	template <>
	const std::string name<quad_2_elem>::value = "Quad 2 elem";
	template <>
	const std::string name<quad_28_elem>::value = "Quad 28 elem";
}


