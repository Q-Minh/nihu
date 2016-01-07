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
 * \file domain.hpp
 * \ingroup funcspace
 * \brief declaration of CRTP base class ::domain_base
 */

#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include "global_definitions.hpp"
#include "space.hpp"

namespace domain_traits
{
	/** \brief assigns a coordinate space to the domain */
    template <class Derived>
    struct space_type;

	/** \brief defines the number of domain corners */
    template <class Derived>
    struct num_corners;

	/** \brief Defines the domain's size (volume) */
	template <class Derived>
	struct volume;

	/** \brief Assigns an id to the domain */
	template <class Derived>
	struct id
	{
		enum {
			value = space_type<Derived>::type::dimension * 10 + num_corners<Derived>::value
		};
	};

	/** \brief Assigns a textual id the domain */
	template <class Derived>
	struct name
	{
		static std::string const value;
	};
}

/** \brief Polygonal subset of the \f$ \xi \f$ space. All elements are defined on a domain. */
template <class Derived>
class domain_base
{
public:
	/** \brief self-returning */
	typedef Derived type;
	/** \brief the space type as nested typedef */
	typedef typename domain_traits::space_type<Derived>::type space_t;
	/** \brief number of domain corners */
	enum { num_corners = domain_traits::num_corners<Derived>::value };
	/** \brief domain id as nested enum */
    enum { id = domain_traits::id<Derived>::value };
    /** \brief space dimensions */
    enum { dimension = space_t::dimension };
	/** \brief coordinate scalar type */
	typedef typename space_t::scalar_t scalar_t;
	/** \brief coordinate vector type */
	typedef typename space_t::location_t xi_t;
	/** \brief type of corners array */
	typedef xi_t corners_t[num_corners];

	/** \brief return domain center */
	static xi_t const &get_center(void)
	{
      	return Derived::get_center_impl();
	}

	/** \brief return begin address of domain corners' array */
	static corners_t const &get_corners(void)
	{
		return Derived::get_corners_impl();
	}

	/** \brief return specified corner of domain */
	static xi_t const &get_corner(unsigned idx)
	{
		return get_corners()[idx];
	}

	/** \brief return domain volume */
	static constexpr scalar_t get_volume(void)
	{
        return domain_traits::volume<Derived>::value;
	}

	/** \brief return domain name */
	static std::string const &get_name(void)
	{
        return domain_traits::name<Derived>::value;
	}

	/** \brief return domain id */
	static constexpr unsigned get_id(void)
	{
        return id;
	}
};

#endif

