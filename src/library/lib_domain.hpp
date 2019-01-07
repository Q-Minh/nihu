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

/** \file lib_domain.hpp
 * \brief implementation of library domains
 */

#ifndef LIB_DOMAIN_HPP_INCLUDED
#define LIB_DOMAIN_HPP_INCLUDED

#include "../core/domain.hpp"

namespace NiHu
{

/** \brief a 1D line domain */
class line_domain;

namespace domain_traits
{
template <>
struct space_type<line_domain> : space_1d<> {};

template <>
struct volume<line_domain> { static constexpr double value = 2.0; };

template <>
struct num_corners<line_domain> { enum { value = 2 }; };

template <>
struct num_edges<line_domain> { enum { value = 1 }; };
}

class line_domain :
	public domain_base<line_domain>
{
public:
	/** \brief return domain corners */
	static corners_t const &get_corners_impl(void) { return m_corners; }
	/** \brief return domain edges */
	static edges_t const &get_edges_impl(void) { return m_edges; }
	/** \brief return domain center */
	static xi_t const &get_center_impl(void) { return m_center; }

private:
	static corners_t const m_corners;
	static edges_t const m_edges;
	static xi_t const m_center;
};

/** \brief a 2D triangle domain */
class tria_domain;

namespace domain_traits
{
template <>
struct space_type<tria_domain> : space_2d<> {};

template <>
struct volume<tria_domain> { static constexpr double value = 0.5; };

template <>
struct num_corners<tria_domain> { enum { value = 3 }; };

template <> 
struct num_edges<tria_domain> { enum { value = 3 }; };
}

class tria_domain :
	public domain_base<tria_domain>
{
public:
	/** \brief return domain corners */
	static corners_t const &get_corners_impl(void) { return m_corners; }
	/** \brief return domain edges */
	static edges_t const &get_edges_impl(void) { return m_edges; }
	/** \brief return domain center */
	static xi_t const &get_center_impl(void) { return m_center; }

private:
	static corners_t const m_corners;
	static edges_t const m_edges;
	static xi_t const m_center;
};

/** \brief a 2D quad domain */
class quad_domain;

namespace domain_traits
{
template <>
struct space_type<quad_domain> : space_2d<> {};

template <> 
struct volume<quad_domain> { static constexpr double value = 4.0; };

template <> 
struct num_corners<quad_domain> { enum { value = 4 }; };

template <>
struct num_edges<quad_domain> { enum { value = 4 }; };
}

class quad_domain :
	public domain_base<quad_domain>
{
public:
	/** \brief return domain corners */
	static corners_t const &get_corners_impl(void) { return m_corners; }
	/** \brief return domain edges */
	static edges_t const &get_edges_impl(void) { return m_edges; }
	/** \brief return domain center */
	static xi_t const &get_center_impl(void) { return m_center; }

private:
	static corners_t const m_corners;
	static edges_t const m_edges;
	static xi_t const m_center;
};

/** \brief a 3D brick domain */
class brick_domain;

namespace domain_traits
{
template <>
struct space_type<brick_domain> : space_3d<> {};

template <>
struct volume<brick_domain> { static constexpr double value = 8.0; };

template <>
struct num_corners<brick_domain> { enum { value = 8 }; };

template <>
struct num_edges<brick_domain> { enum { value = 12 }; };
}

class brick_domain :
	public domain_base<brick_domain>
{
public:
	/** \brief return domain corners */
    static corners_t const &get_corners_impl(void) { return m_corners; }
	/** \brief return domain edges */
    static edges_t const &get_edges_impl(void) { return m_edges; }
	/** \brief return domain center */
    static xi_t const &get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static edges_t const m_edges;
    static xi_t const m_center;
};

} // end of namespace NiHu

#endif // LIB_DOMAIN_HPP_INCLUDED

