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

#ifndef LINE_2_SHAPE_SET_HPP_INCLUDED
#define LINE_2_SHAPE_SET_HPP_INCLUDED

#include "../core/shapeset.hpp"

// Forward declaration
class line_2_shape_set;

namespace shape_set_traits
{
	template <>
	struct domain<line_2_shape_set>
	{
		typedef line_domain type;
	};

	template <>
	struct num_nodes<line_2_shape_set>
	{
		enum { value = 3 };
	};

	template <>
	struct polynomial_order<line_2_shape_set>
	{
		enum { value = 2 };
	};

	template <>
	struct jacobian_order<line_2_shape_set>
	{
		enum { value = 1 };
	};

	template <unsigned Order>
	struct shape_complexity<line_2_shape_set, Order>
	{
		typedef matrix_function_complexity::general type;
	};
}

/** \brief quadratic 3-noded line shape function set */
class line_2_shape_set : public shape_set_base<line_2_shape_set>
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
		int ret = m_domain_indices[idx];
		if (ret < 0)
			throw std::out_of_range("line_2_shape_set domain corner nodes overindexed");
		return ret;
	}

protected:
	/** \brief the corner nodes of the shape set */
	static xi_t const m_corners[num_nodes];

	/** \brief the array of domain corner indices */
	static int const m_domain_indices[num_nodes];
};

line_2_shape_set::xi_t
const line_2_shape_set::m_corners[line_2_shape_set::num_nodes] = {
	line_2_shape_set::xi_t::Constant(-1.0),
	line_2_shape_set::xi_t::Constant(0.0),
	line_2_shape_set::xi_t::Constant(1.0),
};

int const line_2_shape_set::m_domain_indices[line_2_shape_set::num_nodes] = { 0, -1, 1 };


/**
* \brief quadratic 3-noded line shape functions
* \param [in] _xi the domain variable
* \return the shape function vector
*/
template<>
class shape_function<line_2_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<line_2_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<line_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			-xi[0]*(1.0 - xi[0]) / 2.0,
			1.0 - xi[0] * xi[0],
			xi[0] * (1.0 + xi[0]) / 2.0
			).finished();
	}
};

/**
* \brief quadratic 3-noded line shape function derivatives
* \param [in] _xi the domain variable
* \return the shape function gradient matrix
*/
template<>
class shape_function<line_2_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<line_2_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<line_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			xi[0] - 0.5,
			-2.0*xi[0],
			xi[0] + 0.5
			).finished();
	}
};

/**
* \brief quadratic 3-noded line shape function second derivatives
* \return the shape function second derivative matrix
*/
template<>
class shape_function<line_2_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<line_2_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<line_2_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return shape_t(1.0, -2.0, 1.0);
	}
};

#endif // LINE_2_SHAPE_SET_HPP_INCLUDED
