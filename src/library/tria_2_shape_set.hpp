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

#ifndef TRIA_2_SHAPE_SET_HPP_INCLUDED
#define TRIA_2_SHAPE_SET_HPP_INCLUDED

#include "../core/shapeset.hpp"

// Forward declaration
class tria_2_shape_set;

namespace shape_set_traits
{
	template <>
	struct domain<tria_2_shape_set> : tria_domain {};

	template <>
	struct num_nodes<tria_2_shape_set>
	{
		enum { value = 6 };
	};

	template <>
	struct polynomial_order<tria_2_shape_set>
	{
		enum { value = 2 };
	};

	template <>
	struct jacobian_order<tria_2_shape_set>
	{
		enum { value = 2 };
	};

	template <unsigned Order>
	struct shape_complexity<tria_2_shape_set, Order>
	{
		typedef matrix_function_complexity::general type;
	};
}


/**
* \brief quadratic 6-noded tria shape function set
*/
class tria_2_shape_set : public shape_set_base<tria_2_shape_set>
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
			throw std::out_of_range("tria_2_shape_set domain corner nodes overindexed");
		return ret;
	}


protected:
	/** \brief the corner nodes of the shape set */
	static xi_t const m_corners[num_nodes];
	/** \brief the domain's corner indices assigned to the shape set nodes */
	static int const m_domain_corners[num_nodes];
};

tria_2_shape_set::xi_t
const tria_2_shape_set::m_corners[tria_2_shape_set::num_nodes] = {
	tria_2_shape_set::xi_t(0.0, 0.0),
	tria_2_shape_set::xi_t(0.5, 0.0),
	tria_2_shape_set::xi_t(1.0, 0.0),
	tria_2_shape_set::xi_t(0.5, 0.5),
	tria_2_shape_set::xi_t(0.0, 1.0),
	tria_2_shape_set::xi_t(0.0, 0.5)
};

int const tria_2_shape_set::m_domain_corners[tria_2_shape_set::num_nodes] =
{ 0, -1, 1, -1, 2, -1 };

/**
* \brief quadratic 6-noded tria shape functions
* \param [in] _xi the domain variable
* \return the shape function vector
*/
template<>
class shape_function<tria_2_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<tria_2_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<tria_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1];
		return (shape_t() <<
			(eta + xi - 1.0)*(2.0*eta + 2.0*xi - 1.0),
			-4.0*xi*(eta + xi - 1),
			xi*(2.0*xi - 1.0),
			4.0*eta*xi,
			eta*(2.0*eta - 1.0),
			-4.0*eta*(eta + xi - 1.0)
			).finished();
	}
};

/**
* \brief quadratic 6-noded tria shape function derivatives
* \param [in] _xi the domain variable
* \return the shape function gradient matrix
*/
template<>
class shape_function<tria_2_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<tria_2_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<tria_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1];
		return (shape_t() <<
			4.0*eta + 4 * xi - 3, 4.0*eta + 4 * xi - 3.0,
			4.0 - 8 * xi - 4 * eta, -4.0*xi,
			4.0*xi - 1.0, 0.0,
			4.0*eta, 4.0*xi,
			0.0, 4.0*eta - 1.0,
			-4.0*eta, 4.0 - 4.0*xi - 8.0*eta
			).finished();
	}
};

/**
* \brief quadratic 6-noded tria shape function second derivatives
* \return the shape function second derivative matrix
*/
template<>
class shape_function<tria_2_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<tria_2_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<tria_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		return (shape_t() <<
			4.0, 4.0, 4.0,
			-8.0, -4.0, 0.0,
			4.0, 0.0, 0.0,
			0.0, 4.0, 0.0,
			0.0, 0.0, 4.0,
			0.0, -4.0, -8.0
			).finished();
	}
};

#endif // TRIA_2_SHAPE_SET_HPP_INCLUDED
