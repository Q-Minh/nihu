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
 * \file location_normal.hpp
 * \ingroup library
 * \brief implementation of location and normal kernel inputs
 */

#ifndef LOCATION_NORMAL_HPP_INCLUDED
#define LOCATION_NORMAL_HPP_INCLUDED

#include "../util/brick.hpp"
#include "../core/element.hpp"

/** \brief a class representing a simple location brick
 * \tparam Space the coordinate space
 */
template <class Space>
struct location
{
	/**
	 * \brief the brick template class
	 * \tparam wall the wall type to glue the brick on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief template parameter as nested type */
		typedef Space space_t;
		/** \brief the location type */
		typedef typename space_t::location_t x_t;
		/** \brief the scalar type */
		typedef typename space_t::scalar_t scalar_t;

		/** \brief constructor
		 * \tparam elem_t the element type
		 * \param [in] elem the element instance
		 * \param [in] xi the reference domain variable */
		template <class elem_t>
		brick(elem_t const &elem, typename elem_t::xi_t const &xi) :
			wall(elem, xi),
			m_x(elem.get_x(xi))
		{
		}

		/** \brief return the location */
		x_t const &get_x(void) const
		{
			return m_x;
		}

	private:
		x_t m_x;
	};
};


/** \brief a class representing a normal + Jacobian brick
 * \tparam Space the coordinate space
 */
template <class Space>
struct normal_jacobian
{
	/**
	 * \brief the brick template class
	 * \tparam wall the wall type to glue the brick on
	 */
	template <class wall>
	struct brick : public wall
	{
	public:
		/** \brief template parameter as nested type */
		typedef Space space_t;
		/** \brief the location type */
		typedef typename space_t::location_t x_t;
		/** \brief the scalar type */
		typedef typename space_t::scalar_t scalar_t;

		/** \brief constructor
		 * \tparam elem_t the element type
		 * \param [in] elem the element instance
		 * \param [in] xi the reference domain variable */
		template <class elem_t>
		brick(elem_t const &elem, typename elem_t::xi_t const &xi) :
			wall(elem, xi),
			m_norm(elem.get_normal(xi)),
			m_jac(m_norm.norm())
		{
			m_norm /= m_jac;
		}

		/** \brief return the normal */
		x_t const &get_unit_normal(void) const
		{
			return m_norm;
		}

		/** \brief return the Jacobian */
		scalar_t const &get_jacobian(void) const
		{
			return m_jac;
		}

	private:
		x_t m_norm;
		scalar_t m_jac;
	};
};

/** \brief a class representing a Jacobian brick used for volume elements
 * \tparam Space the coordinate space
 */
template <class Space>
struct volume_jacobian
{
	/**
	 * \brief the brick template class
	 * \tparam wall the wall type to glue the brick on
	 */
	template <class wall>
	struct brick : public wall
	{
	public:
		/** \brief template parameter as nested type */
		typedef Space space_t;
		/** \brief the scalar type */
		typedef typename space_t::scalar_t scalar_t;

		/** \brief constructor
		 * \tparam elem_t the element type
		 * \param [in] elem the element instance
		 * \param [in] xi the reference domain variable */
		template <class elem_t>
		brick(elem_t const &elem, typename elem_t::xi_t const &xi) :
			wall(elem, xi),
			m_jac(elem.get_dx(xi).determinant())
		{
		}

		/** \brief return the Jacobian */
		scalar_t const &get_jacobian(void) const
		{
			return m_jac;
		}

	private:
		scalar_t m_jac;
	};
};

typedef build<location<space_1d<> > >::type location_input_1d;
typedef build<location<space_2d<> > >::type location_input_2d;
typedef build<location<space_3d<> > >::type location_input_3d;
typedef build<location<space_2d<> >, normal_jacobian<space_2d<> > >::type location_normal_input_2d;
typedef build<location<space_3d<> >, normal_jacobian<space_3d<> > >::type location_normal_input_3d;

#endif // LOCATION_NORMAL_HPP_INCLUDED

