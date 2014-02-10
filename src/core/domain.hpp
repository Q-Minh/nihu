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
	/** assigns a coordinate space to the domain */
    template <class Derived>
    struct space_type;

	/** assigns the number of domain corners */
    template <class Derived>
    struct num_corners;

	template <class Derived>
	struct volume;

	/** \brief Assign an id to a domain */
	template <class Derived>
	struct id
	{
		enum { value = space_type<Derived>::type::dimension * 10 + num_corners<Derived>::value };
	};

	template <class Derived>
	struct name
	{
		static std::string const value;
	};
}

/**
 * \brief a subset of the \f$\xi\f$ space. All elements are defined on a domain.
 * \tparam Space the coordinate space the domain is defined on
 * \tparam NumCorners the number of corners of the domain
 */
template <class Derived>
class domain_base
{
public:
	typedef Derived type;

	typedef typename domain_traits::space_type<Derived>::type space_t;
	enum
	{
        num_corners = domain_traits::num_corners<Derived>::value,
        id = domain_traits::id<Derived>::value,
        dimension = space_t::dimension
	};

	typedef typename space_t::scalar_t scalar_t;
	typedef typename space_t::location_t xi_t;
	typedef xi_t corners_t[num_corners];

	static xi_t const &get_center(void)
	{
      	static xi_t const center = Derived::get_center_impl();
      	return center;
	}

	static xi_t const *get_corners(void)
	{
		return Derived::get_corners_impl();
	}

	static xi_t const &get_corner(unsigned idx)
	{
		return get_corners()[idx];
	}

	static constexpr scalar_t get_volume(void)
	{
        return domain_traits::volume<Derived>::value;
	}

	static std::string const &get_name(void)
	{
        return domain_traits::name<Derived>::value;
	}
};

#endif

