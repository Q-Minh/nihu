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
 * \brief Definition of various shape function sets
 */

#ifndef SHAPESET_HPP_INCLUDED
#define SHAPESET_HPP_INCLUDED

#include <iostream>
#include <type_traits>
#include <stdexcept>

#include "domain.hpp"
#include "../util/store_or_on_the_fly.hpp"

/** \brief shape function derivative indices */
namespace shape_index
{
	enum {
		dXI = 0,	/**< \brief index of xi in shape derivative matrix */
		dETA = 1,	/**< \brief index of eta in shape derivative matrix */
		dXIXI = 0,	/**< \brief index of xi_xi in shape second derivative matrix */
		dXIETA = 1,	/**< \brief index of xi_eta in shape derivative matrix */
		dETAXI = 1,	/**< \brief index of eta_xi in shape second derivative matrix */
		dETAETA = 2	/**< \brief index of eta-eta in shape second derivative matrix */
	};
}

namespace matrix_function_complexity
{
	struct zero {};
	struct constant {};
	struct general {};
}

template <unsigned Order, unsigned Dim>
struct numnum;

template <unsigned Dim>
struct numnum<0, Dim> { enum { value = 1 }; };

template <unsigned Dim>
struct numnum<1, Dim> { enum { value = Dim }; };

template <unsigned Dim>
struct numnum<2, Dim> { enum { value = Dim * (Dim + 1) / 2 }; };


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

	template <class Derived, unsigned Order>
	struct shape_complexity;

	template <class Derived, unsigned Order>
	struct shape_value_type
	{
		typedef Eigen::Matrix<
			typename domain<Derived>::type::scalar_t,
			num_nodes<Derived>::value,
			numnum<Order, domain<Derived>::type::dimension>::value
		> type;
	};

	template <class Derived, unsigned Order>
	struct shape_return_type : std::conditional<
		std::is_same<
		typename shape_complexity<Derived, Order>::type,
		matrix_function_complexity::general
		>::value,
		typename shape_value_type<Derived, Order>::type,
		typename shape_value_type<Derived, Order>::type const &
	> {};
}

template <class Derived, unsigned Order>
class shape_function;

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
	static typename shape_set_traits::shape_return_type<Derived, Order>::type eval_shape(xi_t const &xi)
	{
		return store_or_on_the_fly<
			std::is_same<
			typename shape_set_traits::shape_complexity<Derived, Order>::type,
			matrix_function_complexity::general
			>::value,
			shape_function<Derived, Order>,
			typename shape_t<Order>::type,
			xi_t
		>::eval(xi);
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
}



/**
* \brief quadratic 3-noded line shape function set
*/
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



// Forward declaration
class tria_2_shape_set;

namespace shape_set_traits
{
	template <>
	struct domain<tria_2_shape_set>
	{
		typedef tria_domain type;
	};

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


// Forward declaration
class quad_2_shape_set;

namespace shape_set_traits
{
	template <>
	struct domain<quad_2_shape_set>
	{
		typedef quad_domain type;
	};

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




#endif // SHAPESET_HPP_INCLUDED

