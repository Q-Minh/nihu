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

/** \file lib_shape.hpp
 * \brief definition of shape function sets
 */

#ifndef LIB_SHAPE_HPP_INCLUDED
#define LIB_SHAPE_HPP_INCLUDED

#include "../core/shapeset.hpp"
#include "lib_domain.hpp"

namespace NiHu
{

/** \brief a constant line shape set */
typedef constant_shape_set<line_domain> line_0_shape_set;
/** \brief a constant triangle shape set */
typedef constant_shape_set<tria_domain> tria_0_shape_set;
/** \brief a constant quadrangle shape set */
typedef constant_shape_set<quad_domain> quad_0_shape_set;
/** \brief a constant brick shape set */
typedef constant_shape_set<brick_domain> brick_0_shape_set;


/** linear line shape set */
typedef isoparam_shape_set<line_domain> line_1_shape_set;

namespace shape_set_traits
{
	template <>
	struct jacobian_order<line_1_shape_set>
	{
		enum { value = 0 };
	};

	template <>
	struct shape_complexity<line_1_shape_set, 0> : matrix_function_complexity::general {};

	template <>
	struct shape_complexity<line_1_shape_set, 1> : matrix_function_complexity::constant {};

	template <>
	struct shape_complexity<line_1_shape_set, 2> : matrix_function_complexity::zero {};
}

/**
 * \brief linear 2-noded line shape functions
 * \param [in] xi domain location vector
 * \return shape function vector set
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (1-\xi)/2 \\ L_2(\xi) = (1+\xi)/2 \f$
 */
template<>
class shape_function<line_1_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<line_1_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<line_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			1.0 - xi[0],
			1.0 + xi[0]
			).finished() / 2.0;
	}
};

/**
* \brief linear 2-noded line shape function derivative matrix
* \return shape function derivative matrix
* \details The shape functions are
*
* \f$L'_1(\xi) = -1/2 \\ L'_2(\xi) = 1/2 \f$
*/
template<>
class shape_function<line_1_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<line_1_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<line_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return (shape_t() <<
			-0.5,
			+0.5
			).finished();
	}
};


/**
* \brief linear 2-noded line shape function second derivative matrix
* \return shape function second derivative matrix
* \details The shape function derivatives are
*
* \f$L''_1(\xi) = 0 \f$
*/
template<>
class shape_function<line_1_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<line_1_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<line_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return shape_t::Zero();
	}
};


/** linear triangle shape set */
typedef isoparam_shape_set<tria_domain> tria_1_shape_set;

namespace shape_set_traits
{
	template <>
	struct jacobian_order<tria_1_shape_set>
	{
		enum { value = 0 };
	};

	template <>
	struct shape_complexity<tria_1_shape_set, 0> : matrix_function_complexity::general {};

	template <>
	struct shape_complexity<tria_1_shape_set, 1> : matrix_function_complexity::constant {};

	template <>
	struct shape_complexity<tria_1_shape_set, 2> : matrix_function_complexity::zero {};
}


/**
* \brief linear 3-noded triangle shape functions
* \param [in] xi the domain variable vector
* \return the shape function vector
* \details The shape functions are
*
* \f$L_1(\xi, \eta) = 1-\xi-\eta \\ L_2(\xi, \eta) = \xi \\ L_3(\xi, \eta) = \eta \f$
*/
template<>
class shape_function<tria_1_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<tria_1_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<tria_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			1.0 - xi[0] - xi[1],
			xi[0],
			xi[1]
			).finished();
	}
};

/**
* \brief linear 3-noded tria elem shape function derivative matrix
* \return sape function derivative matrix
*/
template<>
class shape_function<tria_1_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<tria_1_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<tria_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return (shape_t() <<
			-1.0, -1.0,
			+1.0, 0.0,
			0.0, +1.0
			).finished();
	}
};


/**
* \brief linear 3-noded tria elem shape function second derivative matrix
* \return sape function second derivative matrix
*/
template<>
class shape_function<tria_1_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<tria_1_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<tria_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return shape_t::Zero();
	}
};



/** linear quad shape set */
typedef isoparam_shape_set<quad_domain> quad_1_shape_set;

namespace shape_set_traits
{
	template <>
	struct jacobian_order<quad_1_shape_set>
	{
		enum { value = 1 };
	};

	template <>
	struct shape_complexity<quad_1_shape_set, 0> : matrix_function_complexity::general {};

	template <>
	struct shape_complexity<quad_1_shape_set, 1> : matrix_function_complexity::general {};

	template <>
	struct shape_complexity<quad_1_shape_set, 2> : matrix_function_complexity::constant {};
}


/**
* \brief linear 4-noded general quadrilateral shape functions
* \param [in] xi the domain variable vector
* \return the shape function vector
* \details The shape functions are
*
* \f$L_1(\xi, \eta) = (1-\xi)(1-\eta)/4\\
L_2(\xi, \eta) = (1+\xi)(1-\eta)/4\\
L_3(\xi, \eta) = (1+\xi)(1+\eta)/4\\
L_4(\xi, \eta) = (1-\xi)(1+\eta)/4\f$
*/
template<>
class shape_function<quad_1_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<quad_1_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<quad_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			(1.0 - xi[0])*(1.0 - xi[1]),
			(1.0 + xi[0])*(1.0 - xi[1]),
			(1.0 + xi[0])*(1.0 + xi[1]),
			(1.0 - xi[0])*(1.0 + xi[1])
			).finished() / 4.0;
	}
};

/**
* \brief linear 4-noded general quadrilater shape function derivative matrix
* \return shape function gradient matrix
*/
template<>
class shape_function<quad_1_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<quad_1_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<quad_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			-(1.0 - xi[1]), -(1.0 - xi[0]),
			+(1.0 - xi[1]), -(1.0 + xi[0]),
			+(1.0 + xi[1]), +(1.0 + xi[0]),
			-(1.0 + xi[1]), +(1.0 - xi[0])
			).finished() / 4.0;
	}
};

/**
 * \brief linear 4-noded general quadrilater shape function second derivative matrix
 * \return shape function second derivative matrix
 */
template<>
class shape_function<quad_1_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<quad_1_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<quad_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return (shape_t() <<
			0.0, +.25, 0.0,
			0.0, -.25, 0.0,
			0.0, +.25, 0.0,
			0.0, -.25, 0.0
			).finished();
	}
};


/** linear brick shape set */
typedef isoparam_shape_set<brick_domain> brick_1_shape_set;

namespace shape_set_traits
{
	template <>
	struct jacobian_order<brick_1_shape_set>
	{
		enum { value = 1 };
	};

	template <unsigned Order>
	struct shape_complexity<brick_1_shape_set, Order> : matrix_function_complexity::general {};
}


/**
 * \brief linear 8-noded general brick shape functions
 * \param [in] xi the domain variable vector
 * \return the shape function vector
 * \details The shape functions are
 *
 * \f$L_1(\xi, \eta) = (1-\xi)(1-\eta)(1-\zeta)/8\\
 L_2(\xi, \eta) = (1+\xi)(1-\eta)(1-\zeta)/8\\
 L_3(\xi, \eta) = (1+\xi)(1+\eta)(1-\zeta)/8\\
 L_4(\xi, \eta) = (1-\xi)(1+\eta)(1-\zeta)/8\\
 L_2(\xi, \eta) = (1-\xi)(1-\eta)(1+\zeta)/8\\
 L_2(\xi, \eta) = (1+\xi)(1-\eta)(1+\zeta)/8\\
 L_3(\xi, \eta) = (1+\xi)(1+\eta)(1+\zeta)/8\\
 L_4(\xi, \eta) = (1-\xi)(1+\eta)(1+\zeta)/8\f$
 */
template<>
class shape_function<brick_1_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<brick_1_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<brick_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			(1.0 - xi[0])*(1.0 - xi[1])*(1.0 - xi[2]),
			(1.0 + xi[0])*(1.0 - xi[1])*(1.0 - xi[2]),
			(1.0 + xi[0])*(1.0 + xi[1])*(1.0 - xi[2]),
			(1.0 - xi[0])*(1.0 + xi[1])*(1.0 - xi[2]),
			(1.0 - xi[0])*(1.0 - xi[1])*(1.0 + xi[2]),
			(1.0 + xi[0])*(1.0 - xi[1])*(1.0 + xi[2]),
			(1.0 + xi[0])*(1.0 + xi[1])*(1.0 + xi[2]),
			(1.0 - xi[0])*(1.0 + xi[1])*(1.0 + xi[2])
			).finished() / 8.0;
	}
};

/**
 * \brief linear 8-noded general brick shape function derivative matrix
 * \param [in] xi the domain variable vector
 * \return shape function gradient matrix
 */
template<>
class shape_function<brick_1_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<brick_1_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<brick_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			(-1.0)*(1.0 - xi[1])*(1.0 - xi[2]), (1.0 - xi[0])*(-1.0)*(1.0 - xi[2]), (1.0 - xi[0])*(1.0 - xi[1])*(-1.0),
			(+1.0)*(1.0 - xi[1])*(1.0 - xi[2]), (1.0 + xi[0])*(-1.0)*(1.0 - xi[2]), (1.0 + xi[0])*(1.0 - xi[1])*(-1.0),
			(+1.0)*(1.0 + xi[1])*(1.0 - xi[2]), (1.0 + xi[0])*(+1.0)*(1.0 - xi[2]), (1.0 + xi[0])*(1.0 + xi[1])*(-1.0),
			(-1.0)*(1.0 + xi[1])*(1.0 - xi[2]), (1.0 - xi[0])*(+1.0)*(1.0 - xi[2]), (1.0 - xi[0])*(1.0 + xi[1])*(-1.0),
			(-1.0)*(1.0 - xi[1])*(1.0 + xi[2]), (1.0 - xi[0])*(-1.0)*(1.0 + xi[2]), (1.0 - xi[0])*(1.0 - xi[1])*(+1.0),
			(+1.0)*(1.0 - xi[1])*(1.0 + xi[2]), (1.0 + xi[0])*(-1.0)*(1.0 + xi[2]), (1.0 + xi[0])*(1.0 - xi[1])*(+1.0),
			(+1.0)*(1.0 + xi[1])*(1.0 + xi[2]), (1.0 + xi[0])*(+1.0)*(1.0 + xi[2]), (1.0 + xi[0])*(1.0 + xi[1])*(+1.0),
			(-1.0)*(1.0 + xi[1])*(1.0 + xi[2]), (1.0 - xi[0])*(+1.0)*(1.0 + xi[2]), (1.0 - xi[0])*(1.0 + xi[1])*(+1.0)
			).finished() / 8.0;
	}
};

/**
 * \brief linear 8-noded general brick shape function second derivative matrix
 * \param [in] xi the domain variable vector
 * \return shape function second derivative matrix
 */
template<>
class shape_function<brick_1_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<brick_1_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<brick_1_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			0.0, 1.0 - xi[2], 1.0 - xi[1], 0.0, 1.0 - xi[0], 0.0,
			0.0, xi[2] - 1.0, xi[1] - 1.0, 0.0, xi[0] + 1.0, 0.0,
			0.0, 1.0 - xi[2], -xi[1] - 1.0, 0.0, -xi[0] - 1.0, 0.0,
			0.0, xi[2] - 1.0, xi[1] + 1.0, 0.0, xi[0] - 1.0, 0.0,
			0.0, xi[2] + 1.0, xi[1] - 1.0, 0.0, xi[0] - 1.0, 0.0,
			0.0, -xi[2] - 1.0, 1.0 - xi[1], 0.0, -xi[0] - 1.0, 0.0,
			0.0, xi[2] + 1.0, xi[1] + 1.0, 0.0, xi[0] + 1.0, 0.0,
			0.0, -xi[2] - 1.0, -xi[1] - 1.0, 0.0, 1.0 - xi[0], 0.0
			).finished() / 8.0;
	}
};


// Forward declaration
class parallelogram_shape_set;

namespace shape_set_traits
{
	template <>
	struct domain<parallelogram_shape_set> : quad_domain {};

	template <>
	struct num_nodes<parallelogram_shape_set>
	{
		enum { value = 3 };
	};

	template <>
	struct polynomial_order<parallelogram_shape_set>
	{
		enum { value = 1 };
	};

	template <>
	struct jacobian_order<parallelogram_shape_set>
	{
		enum { value = 0 };
	};

	template <>
	struct shape_complexity<parallelogram_shape_set, 0> : matrix_function_complexity::general {};

	template <>
	struct shape_complexity<parallelogram_shape_set, 1> : matrix_function_complexity::constant {};

	template <>
	struct shape_complexity<parallelogram_shape_set, 2> : matrix_function_complexity::zero {};

	template <>
	struct position_dof_vector<parallelogram_shape_set>
	{
		typedef tmp::vector<dof0, dof0, dof0> type;
	};
}

/**
* \brief linear 3-noded parallelogram shape function set
*/
class parallelogram_shape_set : public shape_set_base<parallelogram_shape_set>
{
public:
	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin_impl(void)
	{
		return domain_t::get_corners();
	}
};

/**
 * \brief linear 3-noded parallelogram shape functions
 * \return the shape function vector
 */
template<>
class shape_function<parallelogram_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<parallelogram_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<parallelogram_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			-xi[0] - xi[1],
			1.0 + xi[0],
			1.0 + xi[1]
			).finished() / 2.0;
	}
};

/**
 * \brief linear 3-noded parallelogram shape function derivatives
 * \return the shape function gradient matrix
 */
template<>
class shape_function<parallelogram_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<parallelogram_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<parallelogram_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return (shape_t() <<
			-1.0, -1.0,
			1.0, 0.0,
			0.0, 1.0
			).finished() / 2.0;
	}
};

/**
 * \brief linear 3-noded parallelogram second shape function derivatives
 * \return the shape function second derivative matrix
 */
template<>
class shape_function<parallelogram_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<parallelogram_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<parallelogram_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return shape_t::Zero();
	}
};


// Forward declaration
class line_2_shape_set;

namespace shape_set_traits
{
	template <>
	struct domain<line_2_shape_set> : line_domain {};

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

	template <>
	struct shape_complexity<line_2_shape_set, 2>
	{
		typedef matrix_function_complexity::constant type;
	};

	template <>
	struct position_dof_vector<line_2_shape_set>
	{
		typedef tmp::vector<dof0, dof1, dof0> type;
	};
}

/** \brief quadratic 3-noded line shape function set */
class line_2_shape_set : public shape_set_base<line_2_shape_set>
{
public:
	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin_impl(void)
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
	static shape_t eval(xi_t const &)
	{
		return shape_t(1.0, -2.0, 1.0);
	}
};


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

	template <>
	struct position_dof_vector<tria_2_shape_set>
	{
		typedef tmp::vector<dof0, dof1, dof0, dof1, dof0, dof1> type;
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
	static xi_t const *corner_begin_impl(void)
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
	static shape_t eval(xi_t const &)
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


// Forward declaration
class quad_2_shape_set;

namespace shape_set_traits
{
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

	template <>
	struct position_dof_vector<quad_2_shape_set>
	{
		typedef tmp::vector<dof0, dof1, dof0, dof1, dof0, dof1, dof0, dof1, dof2> type;
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
	static xi_t const *corner_begin_impl(void)
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

	template <>
	struct position_dof_vector<quad_28_shape_set>
	{
		typedef tmp::vector<dof0, dof1, dof0, dof1, dof0, dof1, dof0, dof1> type;
	};
}

/** \brief quadratic 8-noded quad shape function set */
class quad_28_shape_set : public shape_set_base<quad_28_shape_set>
{
public:
	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin_impl(void)
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

}


#endif // LIB_SHAPE_HPP_INCLUDED

