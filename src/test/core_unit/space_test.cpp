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

#include "core/space.hpp"

// testing if space dimensions are defined correctly
static_assert(space<double, 1>::dimension == 1, "Space1D dimension error");
static_assert(space<double, 2>::dimension == 2, "Space2D dimension error");
static_assert(space<double, 3>::dimension == 3, "Space3D dimension error");
// testing if space scalar is defined correctly
static_assert(std::is_same<space<double, 1>::scalar_t, double>::value, "Space scalar error");
static_assert(std::is_same<space<float, 1>::scalar_t, float>::value, "Space scalar error");
static_assert(std::is_same<space<int, 1>::scalar_t, int>::value, "Space scalar error");
// testing if space shorthands are defined correctly
static_assert(std::is_same<space_1d<float>, space<float, 1> >::value, "Space_1d shorthand error");
static_assert(std::is_same<space_2d<float>, space<float, 2> >::value, "Space_2d shorthand error");
static_assert(std::is_same<space_3d<float>, space<float, 3> >::value, "Space_3d shorthand error");
// testing if default space shorthands are defined correctly
static_assert(std::is_same<space_1d<>, space<double, 1> >::value, "Space_1d default shorthand error");
static_assert(std::is_same<space_2d<>, space<double, 2> >::value, "Space_2d default shorthand error");
static_assert(std::is_same<space_3d<>, space<double, 3> >::value, "Space_3d default shorthand error");

int main() {}
