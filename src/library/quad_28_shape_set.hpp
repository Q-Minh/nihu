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

#ifndef QUAD_28_SHAPE_SET_HPP_INCLUDED
#define QUAD_28_SHAPE_SET_HPP_INCLUDED

#include "../core/shapeset.hpp"


// Forward declaration
class quad_28_shape_set;

namespace shape_set_traits
{
	template <>
	struct domain<quad_28_shape_set> : quad_domain {};

	template <>
	struct num_nodes<quad_28_shape_set>
	{
		enum { value = 8 };
	};

	template <>
	struct polynomial_order<quad_28_shape_set>
	{
		enum { value = 2 };
	};

	template <>
	struct jacobian_order<quad_28_shape_set>
	{
		enum { value = 3 };
	};

	template <unsigned Order>
	struct shape_complexity<quad_28_shape_set, Order>
	{
		typedef matrix_function_complexity::general type;
	};
}

/** \brief quadratic 8-noded quad shape function set */
class quad_28_shape_set : public shape_set_base<quad_28_shape_set>
{
public:
	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return m_corners;
	}

	/** \brief convert a node index to a domain corner index
	 * \param [in] idx the node index
	 * \return the domain corner index.
	 * \details if the domain corner index does not exist, an out_of_range exception is thrown
	 */
	static unsigned node_to_domain_corner(unsigned idx)
	{
		int ret = m_domain_corners[idx];
		if (ret < 0)
			throw std::out_of_range("quad_28_shape_set domain corner nodes overindexed");
		return ret;
	}

protected:
	/** \brief the corner nodes of the shape set */
	static xi_t const m_corners[num_nodes];
	/** \brief the domain corners assigned to the nodal corners */
	static int const m_domain_corners[num_nodes];
};


quad_28_shape_set::xi_t
const quad_28_shape_set::m_corners[quad_28_shape_set::num_nodes] = {
	quad_28_shape_set::xi_t(-1.0, -1.0),
	quad_28_shape_set::xi_t(0.0, -1.0),
	quad_28_shape_set::xi_t(+1.0, -1.0),
	quad_28_shape_set::xi_t(+1.0, 0.0),
	quad_28_shape_set::xi_t(+1.0, +1.0),
	quad_28_shape_set::xi_t(0.0, +1.0),
	quad_28_shape_set::xi_t(-1.0, +1.0),
	quad_28_shape_set::xi_t(-1.0, 0.0)
};

int const quad_28_shape_set::m_domain_corners[quad_28_shape_set::num_nodes] =
{ 0, -1, 1, -1, 2, -1, 3, -1 };


/**
* \brief quadratic 8-noded quad shape functions
* \param [in] _xi the domain variable vector
* \return the shape function vector
*/
template<>
class shape_function<quad_28_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<quad_28_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<quad_28_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1], xi2 = xi*xi, eta2 = eta*eta;
		return (shape_t() <<
			-((xi - 1.0)*(eta - 1.0)*(xi + eta + 1.0)) / 4.0,
			((xi2 - 1.0)*(eta - 1.0)) / 2.0,
			((xi + 1.0)*(eta - 1.0)*(eta - xi + 1.0)) / 4.0,
			-((eta2 - 1.0)*(xi + 1.0)) / 2.0,
			((xi + 1.0)*(eta + 1.0)*(xi + eta - 1.0)) / 4.0,
			-((xi2 - 1.0)*(eta + 1.0)) / 2.0,
			((xi - 1.0)*(eta + 1.0)*(xi - eta + 1.0)) / 4.0,
			((eta2 - 1.0)*(xi - 1.0)) / 2.0
			).finished();
	}
};

/**
* \brief quadratic 8-noded quad shape function derivatives
* \param [in] _xi the domain variable vector
* \return the shape function gradient matrix
*/
template<>
class shape_function<quad_28_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<quad_28_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<quad_28_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto x = _xi[0], y = _xi[1], x2 = x*x, y2 = y*y;
		return (shape_t() <<
			-((2.0*x + y)*(y - 1.0)) / 4.0, -((x + 2 * y)*(x - 1.0)) / 4.0,
			x*(y - 1.0), (x2 - 1.0) / 2.0,
			-((2.0*x - y)*(y - 1.0)) / 4.0, -((x - 2 * y)*(x + 1.0)) / 4.0,
			(1.0 - y2) / 2.0, -y*(x + 1.0),
			((2.0*x + y)*(y + 1.0)) / 4.0, ((x + 2 * y)*(x + 1.0)) / 4.0,
			-x*(y + 1.0), (1.0 - x2) / 2.0,
			((2.0*x - y)*(y + 1.0)) / 4.0, ((x - 2 * y)*(x - 1.0)) / 4.0,
			(y2 - 1.0) / 2.0, y*(x - 1.0)
			).finished();
	}
};

/**
* \brief quadratic 8-noded quad shape function second derivatives
* \param [in] _xi the domain variable vector
* \return the shape function second derivative matrix
*/
template<>
class shape_function<quad_28_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<quad_28_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<quad_28_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi(_xi[0]), eta(_xi[1]);
		return (shape_t() <<
			.5 - eta / 2.0, 1 / 4.0 - xi / 2.0 - eta / 2.0, .5 - xi / 2.0,
			eta - 1.0, xi, 0.0,
			.5 - eta / 2.0, eta / 2.0 - xi / 2.0 - 1 / 4.0, xi / 2.0 + .5,
			0.0, -eta, -xi - 1.0,
			eta / 2.0 + .5, eta / 2.0 + xi / 2.0 + 1 / 4.0, xi / 2.0 + .5,
			-eta - 1.0, -xi, 0.0,
			eta / 2.0 + .5, xi / 2.0 - eta / 2.0 - 1 / 4.0, .5 - xi / 2.0,
			0.0, eta, xi - 1.0
			).finished();
	}
};

#endif
