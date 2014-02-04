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
#include "../util/store_or_on_the_fly.hpp"

/** \brief shape function derivative indices */
namespace shape_derivative_index
{
	enum {
		dXI = 0,	/**< \brief index of xi in 1st derivative matrix */
		dETA = 1,	/**< \brief index of eta in 1st derivative matrix */
		dXIXI = 0,	/**< \brief index of xi_xi in 2nd derivative matrix */
		dXIETA = 1,	/**< \brief index of xi_eta in 2nd derivative matrix */
		dETAXI = 1,	/**< \brief index of eta_xi in 2nd derivative matrix */
		dETAETA = 2	/**< \brief index of eta-eta in 2nd derivative matrix */
	};
}

namespace matrix_function_complexity
{
	struct zero {};
	struct constant {};
	struct general {};
}

template <class Derived, unsigned Order>
class shape_function;

template <unsigned Order, unsigned Dim>
struct num_derivatives;

template <unsigned Dim>
struct num_derivatives<0, Dim> { enum { value = 1 }; };

template <unsigned Dim>
struct num_derivatives<1, Dim> { enum { value = Dim }; };

template <unsigned Dim>
struct num_derivatives<2, Dim> { enum { value = Dim * (Dim + 1) / 2 }; };


/** \brief Traits of shape function sets */
namespace shape_set_traits
{
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
			value = domain_id<typename domain<Derived>::type>::value * 100 +
			num_nodes<Derived>::value
		};
	};

	/** \brief Defines the polynomial order of the shape set */
	template <class Derived>
	struct polynomial_order;

	/** \brief Defines the polynomial order of the shape set's Jacobian */
	template <class Derived>
	struct jacobian_order;

	/** \brief Defines the complexity to determine if the shape functions can be precomputed or not */
	template <class Derived, unsigned Order>
	struct shape_complexity;

	/** \brief Defines the value type of the shape function matrix */
	template <class Derived, unsigned Order>
	struct shape_value_type
	{
		typedef Eigen::Matrix<
			typename domain<Derived>::type::scalar_t,
			num_nodes<Derived>::value,
			num_derivatives<Order, domain<Derived>::type::dimension>::value
		> type;
	};

	/** \brief Defines the factory functor that computes or stores the shape functions */
	template <class Derived, unsigned Order>
	struct factory_functor
	{
		typedef store_or_on_the_fly<
			std::is_same<
			typename shape_complexity<Derived, Order>::type,
			matrix_function_complexity::general
			>::value,
			shape_function<Derived, Order>,
			typename domain<Derived>::type::xi_t
		> type;
	};

	/** \brief Defines the return type of the shape function matrix */
	template <class Derived, unsigned Order>
	struct shape_return_type
	{
		typedef typename factory_functor<Derived, Order>::type::return_type type;
	};
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
	struct shape_t : shape_set_traits::shape_value_type<Derived, Order> {};

	template <unsigned Order>
	static typename shape_set_traits::shape_return_type<Derived, Order>::type
	eval_shape(xi_t const &xi)
	{
		return shape_set_traits::factory_functor<Derived, Order>::type::eval(xi);
	}

	/** \brief return begin iterator of corner nodes */
	static xi_t const *corner_begin(void)
	{
		return Derived::corner_begin();
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
		return *(Derived::corner_begin() + idx);
	}
};

// Forward declaration
template <class Domain>
class constant_shape_set;

namespace shape_set_traits
{
	template <class Domain>
	struct domain<constant_shape_set<Domain> > : Domain{};

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
	struct shape_complexity<constant_shape_set<Domain>, 0>
	{
		typedef matrix_function_complexity::constant type;
	};

	template <class Domain>
	struct shape_complexity<constant_shape_set<Domain>, 1>
	{
		typedef matrix_function_complexity::zero type;
	};

	template <class Domain>
	struct shape_complexity<constant_shape_set<Domain>, 2>
	{
		typedef matrix_function_complexity::zero type;
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
	static xi_t const *corner_begin(void)
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
	static shape_t eval(xi_t const &xi)
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


/** constant line shape set */
typedef constant_shape_set<line_domain> line_0_shape_set;
/** constant triangle shape set */
typedef constant_shape_set<tria_domain> tria_0_shape_set;
/** constant quad shape set */
typedef constant_shape_set<quad_domain> quad_0_shape_set;
/** constant brick shape set */
typedef constant_shape_set<brick_domain> brick_0_shape_set;




// Forward declaration
template <class Domain>
class isoparam_shape_set;

/** \brief traits of shape function sets */
namespace shape_set_traits
{
	template <class Domain>
	struct domain<isoparam_shape_set<Domain> >
	{
		typedef Domain type;
	};

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
	static xi_t const *corner_begin(void)
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
	struct shape_complexity<line_1_shape_set, 0>
	{
		typedef matrix_function_complexity::general type;
	};

	template <>
	struct shape_complexity<line_1_shape_set, 1>
	{
		typedef matrix_function_complexity::constant type;
	};

	template <>
	struct shape_complexity<line_1_shape_set, 2>
	{
		typedef matrix_function_complexity::zero type;
	};
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
	static shape_t eval(xi_t const &xi)
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
	static shape_t eval(xi_t const &xi)
	{
		return shape_t::Zero();
	}
};


/** linear tria shape set */
typedef isoparam_shape_set<tria_domain> tria_1_shape_set;

namespace shape_set_traits
{
	template <>
	struct jacobian_order<tria_1_shape_set>
	{
		enum { value = 0 };
	};

	template <>
	struct shape_complexity<tria_1_shape_set, 0>
	{
		typedef matrix_function_complexity::general type;
	};

	template <>
	struct shape_complexity<tria_1_shape_set, 1>
	{
		typedef matrix_function_complexity::constant type;
	};

	template <>
	struct shape_complexity<tria_1_shape_set, 2>
	{
		typedef matrix_function_complexity::zero type;
	};
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
	static shape_t eval(xi_t const &xi)
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
	static shape_t eval(xi_t const &xi)
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
	struct shape_complexity<quad_1_shape_set, 0>
	{
		typedef matrix_function_complexity::general type;
	};

	template <>
	struct shape_complexity<quad_1_shape_set, 1>
	{
		typedef matrix_function_complexity::general type;
	};

	template <>
	struct shape_complexity<quad_1_shape_set, 2>
	{
		typedef matrix_function_complexity::constant type;
	};
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
	static shape_t eval(xi_t const &xi)
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
	struct shape_complexity<brick_1_shape_set, Order>
	{
		typedef matrix_function_complexity::general type;
	};
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
	struct domain<parallelogram_shape_set>
	{
		typedef quad_domain type;
	};

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
	struct shape_complexity<parallelogram_shape_set, 0>
	{
		typedef matrix_function_complexity::general type;
	};

	template <>
	struct shape_complexity<parallelogram_shape_set, 1>
	{
		typedef matrix_function_complexity::constant type;
	};

	template <>
	struct shape_complexity<parallelogram_shape_set, 2>
	{
		typedef matrix_function_complexity::zero type;
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
	static xi_t const *corner_begin(void)
	{
		return domain_t::get_corners();
	}
};

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
	static shape_t eval(xi_t const &xi)
	{
		return (shape_t() <<
			-.5, -.5,
			+.5, .0,
			.0, +.5
			).finished();
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
	static shape_t eval(xi_t const &xi)
	{
		return shape_t::Zero();
	}
};

#endif // SHAPESET_HPP_INCLUDED
