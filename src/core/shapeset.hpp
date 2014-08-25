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
 * \file shapeset.hpp
 * \ingroup funcspace
 * \brief Definition of shape function sets
 */

#ifndef SHAPESET_HPP_INCLUDED
#define SHAPESET_HPP_INCLUDED

#include <iostream>
#include <type_traits>
#include <stdexcept>

#include "domain.hpp"
#include "../tmp/vector.hpp"
#include "../tmp/algorithm.hpp"
#include "../util/conditional_precompute.hpp"

/** \brief shape function derivative indices
 * These indices are used to index the first and second derivatives of the
 * shape functions.
 */
namespace shape_derivative_index
{
	enum {
		dXI = 0,	/**< \brief index of xi in 1st derivative matrix */
		dETA = 1,	/**< \brief index of eta in 1st derivative matrix */
		dXIXI = 0,	/**< \brief index of xi_xi in 2nd derivative matrix */
		dXIETA = 1,	/**< \brief index of xi_eta in 2nd derivative matrix */
		dETAXI = 1,	/**< \brief index of eta_xi in 2nd derivative matrix */
		dETAETA = 2	/**< \brief index of eta_eta in 2nd derivative matrix */
	};
}

/** \brief position degree of freedom of a point in the intrinsic domain
 * \details Position degree of freedom is the number of independent directions
 * where the intrinsic domain is open.
 * For example, a corner point in a 2D quad domain is of 0DOF,
 * a point on the edge is 1DOF, and an internal point is 2DOF.
 */
template <unsigned d>
struct position_dof : std::integral_constant<unsigned, d> {};
/** \brief shorthand for 0 dof */
typedef position_dof<0> dof0;
/** \brief shorthand for 1 dof */
typedef position_dof<1> dof1;
/** \brief shorthand for 2 dof */
typedef position_dof<2> dof2;

// forward declaration
template <class Derived, unsigned Order>
class shape_function;

/** \brief return the number of partial derivatives in dim dimensions
 * In 1D there is one partial derivative and one second partial derivative.
 * In 2D there are two partial derivatives and four second derivatives, but
 * because of symmetry, the actual number of different second derivatives is 3.
 */
constexpr unsigned num_derivatives(unsigned order, unsigned dim)
{
	return order == 0 ? 1 : (order == 1 ? dim : dim * (dim + 1) / 2);
}

/** \brief Traits of shape function sets */
namespace shape_set_traits
{
	/** \brief The shape set's textual id - used for debug information */
	template <class Derived>
	struct name
	{
		static const std::string value;	/**< \brief the textual name */
	};

	/** \brief Defines the domain where the shape function set is defined */
	template <class Derived>
	struct domain;

	/** \brief Defines the number of shape functions in the set */
	template <class Derived>
	struct num_nodes;

	/** \brief Assigns an id to the shape set */
	template <class Derived>
	struct id
	{
		enum {
			value = domain_traits::id<typename domain<Derived>::type>::value * 100 +
			num_nodes<Derived>::value
		};
	};

	/** \brief Defines the polynomial order of the shape set */
	template <class Derived>
	struct polynomial_order;

	/** \brief Defines the polynomial order of the shape set's Jacobian
	 * \todo this should be moved to element_traits, where the coordinate transform is defined
	 */
	template <class Derived>
	struct jacobian_order;

	/** \brief Defines the complexity to determine if the shape functions can be precomputed or not */
	template <class Derived, unsigned Order>
	struct shape_complexity;

	/** \brief Defines the value type of the shape function matrix (and derivatives) */
	template <class Derived, unsigned Order>
	struct shape_value_type
	{
		typedef Eigen::Matrix<
			typename domain<Derived>::type::scalar_t,
			num_nodes<Derived>::value,
			num_derivatives(Order, domain<Derived>::type::dimension)
		> type;
	};

	/** \brief Defines the factory functor that computes or stores the shape functions */
	template <class Derived, unsigned Order>
	struct factory_functor : conditional_precompute<
		typename shape_complexity<Derived, Order>::type,
		shape_function<Derived, Order>,
		typename domain<Derived>::type::xi_t
	> {};

	/** \brief Defines the return type of the shape function matrix */
	template <class Derived, unsigned Order>
	struct shape_return_type
	{
		typedef typename factory_functor<Derived, Order>::type::return_type type;
	};

	/** \brief defines the nodal degrees of freedoms of the shape functions */
	template <class Derived>
	struct position_dof_vector;
}


/**
 * \brief Shapeset base class for CRTP
 * \tparam Derived The derived shapeset class
 */
template <class Derived>
class shape_set_base
{
public:
	/** \brief Domain */
	typedef typename shape_set_traits::domain<Derived>::type domain_t;

	/** \brief integer constants */
	enum {
		/** number of shape function nodes */
		num_nodes = shape_set_traits::num_nodes<Derived>::value,
		/** \brief the shape set polynomial order */
		polynomial_order = shape_set_traits::polynomial_order<Derived>::value,
		/** \brief the Jacobian polynomial order */
		jacobian_order = shape_set_traits::jacobian_order<Derived>::value,
		/** \brief the shape set id */
		id = shape_set_traits::id<Derived>::value,
	};

	/** \brief scalar type inherited from the domain */
	typedef typename domain_t::scalar_t scalar_t;
	/** type of the local coordinate */
	typedef typename domain_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector */
	template <unsigned Order>
	struct shape_value_type : shape_set_traits::shape_value_type<Derived, Order> {};

	typedef typename shape_value_type<0>::type shape_t;
	typedef typename shape_value_type<1>::type dshape_t;
	typedef typename shape_value_type<2>::type ddshape_t;

	template <unsigned Order>
	static typename shape_set_traits::shape_return_type<Derived, Order>::type
	eval_shape(xi_t const &xi)
	{
		return shape_set_traits::factory_functor<Derived, Order>::type::eval(xi);
	}

	typedef typename shape_set_traits::position_dof_vector<Derived>::type position_dof_vector;

	/** \brief return begin iterator of corner nodes */
	static xi_t const *corner_begin(void)
	{
		return Derived::corner_begin_impl();
	}

	/** \brief return end iterator of corner nodes */
	static xi_t const *corner_end(void)
	{
		return Derived::corner_begin() + num_nodes;
	}

	/** \brief return corner at a given node number
	 * \param [in] idx the node index
	 * \return constant reference to the coordinates
	 */
	static xi_t const &corner_at(size_t idx)
	{
		return *(corner_begin() + idx);
	}
};

// Forward declaration
template <class Domain>
class constant_shape_set;

namespace shape_set_traits
{
	template <class Domain>
	struct domain<constant_shape_set<Domain> > : Domain {};

	template <class Domain>
	struct num_nodes<constant_shape_set<Domain> >
	{
		enum { value = 1 };
	};

	template <class Domain>
	struct polynomial_order<constant_shape_set<Domain> >
	{
		enum { value = 0 };
	};

	template <class Domain>
	struct jacobian_order<constant_shape_set<Domain> >
	{
		enum { value = 0 };
	};

	template <class Domain>
	struct shape_complexity<constant_shape_set<Domain>, 0> : matrix_function_complexity::constant {};

	template <class Domain>
	struct shape_complexity<constant_shape_set<Domain>, 1> : matrix_function_complexity::zero {};

	template <class Domain>
	struct shape_complexity<constant_shape_set<Domain>, 2> : matrix_function_complexity::zero {};

	template <class Domain>
	struct position_dof_vector<constant_shape_set<Domain> >
	{
		typedef tmp::vector<position_dof<Domain::dimension> > type;
	};
}

/**
* \brief Constant interpolation functions
* \tparam Domain the reference domain the shape set is defined above
*/
template <class Domain>
class constant_shape_set : public shape_set_base<constant_shape_set<Domain> >
{
public:
	/** \brief the CRTP base class type */
	typedef shape_set_base<constant_shape_set<Domain> > base_t;
	/** \brief the domain type */
	typedef typename base_t::domain_t domain_t;
	/** \brief type of the domain variable */
	typedef typename base_t::xi_t xi_t;

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin_impl(void)
	{
		return &(domain_t::get_center());
	}
};


template <class Domain>
class shape_function<constant_shape_set<Domain>, 0>
{
	typedef typename shape_set_traits::shape_value_type<constant_shape_set<Domain>, 0>::type shape_t;
	typedef typename shape_set_traits::domain<constant_shape_set<Domain> >::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return shape_t::Ones();
	}
};

template <class Domain>
class shape_function<constant_shape_set<Domain>, 1>
{
	typedef typename shape_set_traits::shape_value_type<constant_shape_set<Domain>, 1>::type shape_t;
	typedef typename shape_set_traits::domain<constant_shape_set<Domain>>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return shape_t::Zero();
	}
};

template <class Domain>
class shape_function<constant_shape_set<Domain>, 2>
{
	typedef typename shape_set_traits::shape_value_type<constant_shape_set<Domain>, 2>::type shape_t;
	typedef typename shape_set_traits::domain<constant_shape_set<Domain>>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &xi)
	{
		return shape_t::Zero();
	}
};


// Forward declaration
template <class Domain>
class isoparam_shape_set;

/** \brief traits of isoparametric shape function sets */
namespace shape_set_traits
{
	template <class Domain>
	struct domain<isoparam_shape_set<Domain> > : Domain {};

	template <class Domain>
	struct num_nodes<isoparam_shape_set<Domain> >
	{
		enum { value = Domain::num_corners };
	};

	template <class Domain>
	struct polynomial_order<isoparam_shape_set<Domain> >
	{
		enum { value = 1 };
	};

	template <class Domain>
	struct position_dof_vector<isoparam_shape_set<Domain> > : tmp::constant_sequence<
		dof0,
		num_nodes<isoparam_shape_set<Domain> >::value,
		tmp::vector<>
	> {};
}

/**
* \brief Isoparametric shape sets
*/
template <class Domain>
class isoparam_shape_set : public shape_set_base<isoparam_shape_set<Domain> >
{
public:
	/** \brief the CRTP base class */
	typedef shape_set_base<isoparam_shape_set<Domain> > base_t;
	/** \brief the domain type */
	typedef typename base_t::domain_t domain_t;

	/** \brief the xi location vector type */
	typedef typename base_t::xi_t xi_t;

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin_impl(void)
	{
		return domain_t::get_corners();
	}

	/** \brief assign domain corner to node
	 * \param [in] idx node index
	 * \return the same domain index
	 */
	static unsigned node_to_domain_corner(unsigned idx)
	{
		return idx;
	}
};

#endif // SHAPESET_HPP_INCLUDED

