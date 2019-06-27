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

#include "../core/element.hpp"

namespace NiHu
{

/** \brief a class representing a simple location
 * \tparam Space the coordinate space
 */
template <class Space>
class location_input
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
	location_input(elem_t const &elem, typename elem_t::xi_t const &xi) :
		m_x(elem.get_x(xi))
	{
	}

	/** \brief constructor from location
	 * \param [in] x the location */
	explicit location_input(x_t const &x)
		: m_x(x)
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


/** \brief a class representing a normal + Jacobian input
 * \tparam Space the coordinate space
 */
template <class Space>
class location_normal_jacobian_input
	: public location_input<Space>
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
	location_normal_jacobian_input(elem_t const &elem, typename elem_t::xi_t const &xi)
		: location_input<Space>(elem, xi)
		, m_jac_vector(elem.get_normal(xi))
		, m_jac(m_jac_vector.norm())
		, m_unit_normal(m_jac_vector / m_jac)
	{
	}

	/** \brief constructor from location and jacobian vector
	 * \param [in] x the location
	 * \param [in] jac the jacobian vector
	 */
	location_normal_jacobian_input(x_t const &x, x_t const &jac_vector)
		: location_input<Space>(x)
		, m_jac_vector(jac_vector)
		, m_jac(m_jac_vector.norm())
		, m_unit_normal(m_jac_vector / m_jac)
	{
	}

	/** \brief return the unit normal vector */
	x_t const &get_unit_normal(void) const
	{
		return m_unit_normal;
	}

	/** \brief return the Jacobian */
	scalar_t const &get_jacobian(void) const
	{
		return m_jac;
	}
	
	/// \brief return the jacobian
	/// \return the jacobian
	x_t const &get_jacobian_vector() const
	{
		return m_jac_vector;
	}

private:
	x_t m_jac_vector;
	scalar_t m_jac;
	x_t m_unit_normal;
};


/** \brief a class representing a Jacobian input used for volume elements
 * \tparam Space the coordinate space
 */
template <class Space>
class location_volume_jacobian_input
	: public location_input<Space>
{
public:
	/** \brief template parameter as nested type */
	typedef Space space_t;
	/** \brief the scalar type */
	typedef typename space_t::scalar_t scalar_t;
	
	typedef typename space_t::location_t location_t;

	/** \brief constructor
		* \tparam elem_t the element type
		* \param [in] elem the element instance
		* \param [in] xi the reference domain variable */
	template <class elem_t>
	location_volume_jacobian_input(elem_t const &elem, typename elem_t::xi_t const &xi) :
		location_input<Space>(elem, xi),
		m_jac(elem.get_dx(xi).determinant())
	{
	}

	/// \brief constructor
	/// \param [in] x the location
	/// \param [in] jac the jacobian
	location_volume_jacobian_input(location_t const &x, scalar_t const &jac) :
		location_input<Space>(x),
		m_jac(jac)
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


typedef location_input<space_1d<> > location_input_1d;
typedef location_input<space_2d<> > location_input_2d;
typedef location_input<space_3d<> > location_input_3d;

typedef location_normal_jacobian_input<space_2d<> > location_normal_input_2d;
typedef location_normal_jacobian_input<space_3d<> > location_normal_input_3d;


} // end of namespace NiHu

#endif // LOCATION_NORMAL_HPP_INCLUDED

