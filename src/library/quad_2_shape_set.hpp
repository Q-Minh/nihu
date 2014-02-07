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

#ifndef QUAD_2_SHAPE_SET_HPP_INCLUDED
#define QUAD_2_SHAPE_SET_HPP_INCLUDED

#include "../core/shapeset.hpp"

// Forward declaration
class quad_2_shape_set;

namespace shape_set_traits
{
	template <>
	const std::string name<quad_2_shape_set>::value = "Quad 2 shape set";

	template <>
	struct domain<quad_2_shape_set> : quad_domain {};

	template <>
	struct num_nodes<quad_2_shape_set>
	{
		enum { value = 9 };
	};

	template <>
	struct polynomial_order<quad_2_shape_set>
	{
		enum { value = 2 };
	};

	template <>
	struct jacobian_order<quad_2_shape_set>
	{
		enum { value = 3 };
	};

	template <unsigned Order>
	struct shape_complexity<quad_2_shape_set, Order>
	{
		typedef matrix_function_complexity::general type;
	};
}


/**
* \brief quadratic 9-noded quad shape function set
*/
class quad_2_shape_set : public shape_set_base<quad_2_shape_set>
{
public:
	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return m_corners;
	}

	/** \brief assign a domain corner to a shapeset node corner
	 * \param [in] idx the indes of the shapeset node
	 * \return index of the domain node
	 * \details the function throws an exception if nonexisting index is searched
	 */
	static unsigned node_to_domain_corner(unsigned idx)
	{
		int ret = m_domain_corners[idx];
		if (ret < 0)
			throw std::out_of_range("quad_2_shape_set domain corner nodes overindexed");
		return ret;
	}

protected:
	/** \brief the corner nodes of the shape set */
	static xi_t const m_corners[num_nodes];
	/** \brief the domain's corner indices assigned to the shape set nodes */
	static int const m_domain_corners[num_nodes];
};


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


/**
* \brief quadratic 9-noded quad shape functions
* \param [in] _xi the domain variable
* \return the shape function vector
*/
template<>
class shape_function<quad_2_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<quad_2_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<quad_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1];
		auto _1mxi = 1 - xi, _1pxi = 1 + xi;
		auto _1meta = 1 - eta, _1peta = 1 + eta;
		return (shape_t() <<
			_1mxi*xi * _1meta*eta / 4.0,
			_1mxi*_1pxi * _1meta*(-eta) / 2.0,
			_1pxi*xi * _1meta*(-eta) / 4.0,
			_1pxi*xi * _1meta*_1peta / 2.0,
			_1pxi*xi * _1peta*eta / 4.0,
			_1mxi*_1pxi * _1peta*eta / 2.0,
			_1mxi*(-xi) * _1peta*eta / 4.0,
			_1mxi*(-xi) * _1meta*_1peta / 2.0,
			_1mxi*_1pxi * _1meta*_1peta
			).finished();
	}
};

/**
* \brief quadratic 9-noded quad shape function derivatives
* \param [in] _xi the domain variable
* \return the shape function gradient matrix
*/
template<>
class shape_function<quad_2_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<quad_2_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<quad_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1], xi2 = xi*xi, eta2 = eta*eta;
		return (shape_t() <<
			eta*(2.0*xi - 1.0)*(eta - 1.0) / 4.0, xi*(2.0*eta - 1.0)*(xi - 1.0) / 4.0,
			-xi*eta*(eta - 1.0), -(xi2 - 1.0)*(2.0*eta - 1.0) / 2.0,
			eta*(2.0*xi + 1.0)*(eta - 1.0) / 4.0, xi*(2.0*eta - 1.0)*(xi + 1.0) / 4.0,
			-(2.0*xi + 1.0)*(eta2 - 1.0) / 2.0, -xi*eta*(xi + 1.0),
			eta*(2.0*xi + 1.0)*(eta + 1.0) / 4.0, xi*(2.0*eta + 1.0)*(xi + 1.0) / 4.0,
			-xi*eta*(eta + 1.0), -(xi2 - 1.0)*(2.0*eta + 1.0) / 2.0,
			eta*(2.0*xi - 1.0)*(eta + 1.0) / 4.0, xi*(2.0*eta + 1.0)*(xi - 1.0) / 4.0,
			-(2.0*xi - 1.0)*(eta2 - 1.0) / 2.0, -xi*eta*(xi - 1.0),
			2.0*xi*(eta2 - 1.0), 2.0*eta*(xi2 - 1.0)
			).finished();
	}
};

/**
 * \brief quadratic 9-noded quad shape function second derivatives
 * \param [in] _xi the domain variable
 * \return the shape function second derivative matrix
 */
template<>
class shape_function<quad_2_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<quad_2_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<quad_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1], xi2 = xi*xi, eta2 = eta*eta;
		return (shape_t() <<
			eta2 / 2.0 - eta / 2.0, eta*xi - xi / 2.0 - eta / 2.0 + 1.0 / 4.0, xi2 / 2.0 - xi / 2.0,
			-eta2 + eta, xi - 2 * eta*xi, 1.0 - xi2,
			eta2 / 2.0 - eta / 2.0, eta / 2.0 - xi / 2.0 + eta*xi - 1.0 / 4.0, xi2 / 2.0 + xi / 2.0,
			1.0 - eta2, -eta - 2 * eta*xi, -xi2 - xi,
			eta2 / 2.0 + eta / 2.0, eta / 2.0 + xi / 2.0 + eta*xi + 1.0 / 4.0, xi2 / 2.0 + xi / 2.0,
			-eta2 - eta, -xi - 2 * eta*xi, 1.0 - xi2,
			eta2 / 2.0 + eta / 2.0, xi / 2.0 - eta / 2.0 + eta*xi - 1.0 / 4.0, xi2 / 2.0 - xi / 2.0,
			1.0 - eta2, eta - 2 * eta*xi, -xi2 + xi,
			2 * eta2 - 2.0, 4 * eta*xi, 2 * xi2 - 2.0
			).finished();
	}
};

#endif // QUAD_2_SHAPE_SET_HPP_INCLUDED
