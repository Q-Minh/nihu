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
 * \file domain.hpp
 * \ingroup funcspace
 * \brief declaration of class ::domain
 */

#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include "global_definitions.hpp"
#include "space.hpp"

namespace domain_traits
{
	/** \brief Assign a name to the domain */
	template <class Domain>
	struct name
	{
		static std::string const value;
	};

	/** \brief Assign an id to a domain */
	template <class Domain>
	struct id
	{
		enum { value = Domain::dimension * 10 + Domain::num_corners };
	};
}

/**
 * \brief a subset of the \f$\xi\f$ space. All elements are defined on a domain.
 * \tparam Space the coordinate space the domain is defined on
 * \tparam NumCorners the number of corners of the domain
 */
template <class Space, unsigned NumCorners>
class domain
{
public:
	/** \brief self-returning */
	typedef domain type;

	/** \brief template parameter as nested type */
	typedef Space space_t;

	/** \brief scalar type inherited from the space */
	typedef typename space_t::scalar_t scalar_t;
	/** \brief dimension inherited from the space */
	static unsigned const dimension = space_t::dimension;

	enum {
		/** \brief number of domain corners */
		num_corners = NumCorners,
		/** \brief the domain id */
		id = domain_traits::id<domain>::value
	};

	/** \brief location vector renamed */
	typedef typename space_t::location_t xi_t;

	/** \brief type of the corners' array */
	typedef const xi_t corners_t[num_corners];

	/**
	 * \brief return the central point of the domain
	 * \return center point
	 */
	static constexpr xi_t const &get_center(void)
	{
		return m_center;
	}

	/**
	 * \brief return pointer to array of corners
	 * \return pointer to array of corners
	 */
	static constexpr xi_t const *get_corners(void)
	{
		return m_corners;
	}

	/**
	 * \brief return reference to a corner point
	 * \return constant reference to a corner point
	 */
	static constexpr xi_t const &get_corner(unsigned idx)
	{
		return get_corners()[idx];
	}

	/**
	 * \brief return domain volume
	 * \return domain volume
	 */
	static constexpr scalar_t const &get_volume(void)
	{
		return m_volume;
	}

protected:
	/** \brief the central point of the domain */
	static xi_t const m_center;
	/** \brief the corner points of the domain */
	static corners_t const m_corners;
	/** \brief the domain's volume */
	static scalar_t const m_volume;
};


typedef domain<space_1d, 2> line_domain;
typedef domain<space_2d, 3> tria_domain;
typedef domain<space_2d, 4> quad_domain;
typedef domain<space_3d, 8> brick_domain;

namespace domain_traits
{
	template <>	std::string const name<line_domain>::value = "1D Line domain";
	template <>	std::string const name<tria_domain>::value = "2D Tria domain";
	template <>	std::string const name<quad_domain>::value = "2D Quad domain";
	template <>	std::string const name<brick_domain>::value = "3D Brick domain";
}

template<>
line_domain::xi_t
	const line_domain::m_center =
	line_domain::xi_t::Zero();

template<>
line_domain::corners_t
	const line_domain::m_corners = {
	line_domain::xi_t::Constant(-1.0),
	line_domain::xi_t::Constant(1.0)
	};

template <>
line_domain::scalar_t
	const line_domain::m_volume = 2.0;

template<>
tria_domain::xi_t
	const tria_domain::m_center =
	tria_domain::xi_t::Constant(1.0/3.0);

template<>
tria_domain::corners_t const tria_domain::m_corners = {
	tria_domain::xi_t(0.0,0.0),
	tria_domain::xi_t(1.0,0.0),
	tria_domain::xi_t(0.0,1.0)
	};

template <>
tria_domain::scalar_t
	const tria_domain::m_volume = 0.5;

template<>
quad_domain::xi_t
	const quad_domain::m_center =
	quad_domain::xi_t::Zero();

template<>
quad_domain::corners_t const quad_domain::m_corners = {
	quad_domain::xi_t(-1.0,-1.0),
	quad_domain::xi_t( 1.0,-1.0),
	quad_domain::xi_t( 1.0, 1.0),
	quad_domain::xi_t(-1.0, 1.0)
	};

template <>
quad_domain::scalar_t
	const quad_domain::m_volume = 4.0;

template<>
brick_domain::corners_t
	const brick_domain::m_corners = {
	brick_domain::xi_t(-1.0,-1.0,-1.0),
	brick_domain::xi_t( 1.0,-1.0,-1.0),
	brick_domain::xi_t( 1.0, 1.0,-1.0),
	brick_domain::xi_t(-1.0, 1.0,-1.0),
	brick_domain::xi_t(-1.0,-1.0, 1.0),
	brick_domain::xi_t( 1.0,-1.0, 1.0),
	brick_domain::xi_t( 1.0, 1.0, 1.0),
	brick_domain::xi_t(-1.0, 1.0, 1.0)
};

template<>
brick_domain::xi_t
	const brick_domain::m_center =
	brick_domain::xi_t::Zero();

template <>
brick_domain::scalar_t
	const brick_domain::m_volume = 8.0;

#endif

