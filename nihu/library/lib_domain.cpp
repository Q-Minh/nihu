// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2019  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2019  Peter Rucz <rucz@hit.bme.hu>
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
 * \file lib_domain.cpp
 * \brief domain library
 * \details
 * This library contains the supported domain types and their constant
 * properties.
 */

#include "lib_domain.hpp"

namespace NiHu
{

// definition of line_domain

namespace domain_traits
{
template <>
const std::string name<line_domain>::value = "Line domain";
}

line_domain::corners_t
const line_domain::m_corners = {
line_domain::xi_t::Constant(-1.),
line_domain::xi_t::Constant(1.)
};

line_domain::edges_t const line_domain::m_edges = { {0, 1} };

line_domain::xi_t const line_domain::m_center = line_domain::xi_t::Zero();


// definition of tria_domain

namespace domain_traits
{
template <>
std::string const name<tria_domain>::value = "Tria domain";
}

tria_domain::corners_t
const tria_domain::m_corners = {
tria_domain::xi_t(0.,0.),
tria_domain::xi_t(1.,0.),
tria_domain::xi_t(0.,1.)
};

tria_domain::edges_t const tria_domain::m_edges = {
	{0, 1},  {1, 2},  {2, 0}
};

tria_domain::xi_t const tria_domain::m_center
= tria_domain::xi_t::Constant(1. / 3.);


// definition of quad_domain

namespace domain_traits
{
template <>
std::string const name<quad_domain>::value = "Quad domain";
}

quad_domain::corners_t
const quad_domain::m_corners = {
quad_domain::xi_t(-1.,-1.),
quad_domain::xi_t(1.,-1.),
quad_domain::xi_t(1., 1.),
quad_domain::xi_t(-1., 1.)
};

quad_domain::edges_t const quad_domain::m_edges = {
	{0, 1},  {1, 2},  {2, 3}, {3, 0}
};

quad_domain::xi_t const quad_domain::m_center = quad_domain::xi_t::Zero();


// definition of brick_domain

namespace domain_traits
{
template <>
std::string const name<brick_domain>::value = "Brick domain";
}

brick_domain::corners_t
	const brick_domain::m_corners = {
	brick_domain::xi_t(-1.,-1.,-1.), brick_domain::xi_t( 1.,-1.,-1.),
	brick_domain::xi_t( 1., 1.,-1.), brick_domain::xi_t(-1., 1.,-1.),
	brick_domain::xi_t(-1.,-1., 1.), brick_domain::xi_t( 1.,-1., 1.),
	brick_domain::xi_t( 1., 1., 1.), brick_domain::xi_t(-1., 1., 1.)
};

brick_domain::edges_t const brick_domain::m_edges = {
	{0, 1},  {1, 2},  {2, 3}, {3, 0},
	{4+0, 4+1},  {4+1, 4+2},  {4+2, 4+3}, {4+3, 4+0},
	{0, 4},  {1, 5},  {2, 6}, {3, 7}
 };

brick_domain::xi_t const brick_domain::m_center = brick_domain::xi_t::Zero();

} // end of namespace NiHu
