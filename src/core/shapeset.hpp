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
#include <stdexcept>

#include "domain.hpp"

/** \brief Traits for shapesets */
template <class Dervied>
struct shape_set_traits;

/** \brief metafunction assigning a shape set id to a shape set
* \brief tparam shape_set the shape set
*/
template <class shape_set>
struct shape_set_id
{
	/** \brief default shape set id computation */
	enum { value =
		domain_id<typename shape_set_traits<shape_set>::domain_t>::value * 100 +
		shape_set_traits<shape_set>::num_nodes
	};
};

/**
* \brief Shapeset base class for CRTP
* \tparam Derived The derived shapeset class
*/
template <class Derived>
class shape_set_base
{
public:
	/** \brief the traits class type */
	typedef shape_set_traits<Derived> traits_t;

	/** \brief Domain */
	typedef typename traits_t::domain_t domain_t;

	/** \brief integer constants */
	enum {
		/** number of shape function nodes */
		num_nodes = traits_t::num_nodes,
		/** \brief the shape set polynomial order */
		polynomial_order = traits_t::polynomial_order,
		/** \brief the Jacobian polynomial order */
		jacobian_order = traits_t::jacobian_order,
		/** \brief the shape set id */
		id = shape_set_id<Derived>::value,
		/** \brief number of second derivatives */
		num_dd = domain_t::dimension * (domain_t::dimension+1) / 2
	};

	/** \brief scalar type inherited from the domain */
	typedef typename domain_t::scalar_t scalar_t;
	/** type of the local coordinate */
	typedef typename domain_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector */
	typedef Eigen::Matrix<scalar_t, num_nodes, 1> shape_t;
	/** \brief type of a \f$\nabla L(\xi)\f$ gradient matrix */
	typedef Eigen::Matrix<scalar_t, num_nodes, domain_t::dimension> dshape_t;
	/** \brief type of the double derivative \f$ L''(\xi)\f$ matrix */
	typedef Eigen::Matrix<scalar_t, num_nodes, num_dd> ddshape_t;

public:
	/** \brief return end iterator to the corner nodes
	* \return end iterator to corner nodes
	*/
	static xi_t const *corner_end(void)
	{
		return Derived::corner_begin() + num_nodes;
	}

	/**
	* \brief return corner at a given node number
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

/** \brief Traits for constant shape sets */
template <class Domain>
struct shape_set_traits<constant_shape_set<Domain> >
{
	/** \brief the domain type */
	typedef Domain domain_t;
	/** \brief number of nodes */
	enum {
		num_nodes = 1,
		polynomial_order = 0,
		jacobian_order = 0
	};
};

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
	/** \brief type of a shape function vector */
	typedef typename base_t::shape_t shape_t;
	/** \brief type of a shape function gradient matrix */
	typedef typename base_t::dshape_t dshape_t;
	/** \brief type of a shape function second derivative matrix */
	typedef typename base_t::ddshape_t ddshape_t;

	/**
	* \brief return constant shape functions
	* \return shape function matrix
	* \details The shape functions are
	*
	* \f$L_1(\xi) = 1 \f$
	*/
	static shape_t const &eval_shape(xi_t const &)
	{
		return m_shape;
	}

	/**
	* \brief Derivatives of constant shape functions
	* \return derivatives of the constant shape functions
	* \details The derivatives are
	*
	* \f$\nabla L_1(\xi) = 0\f$
	*/
	static dshape_t const &eval_dshape(xi_t const &)
	{
		return m_dshape;
	}

	/**
	* \brief Second derivatives of constant shape functions
	* \return Second derivatives of the constant shape functions
	* \details The derivatives are
	*
	* \f$\L''_1(\xi) = 0\f$
	*/
	static ddshape_t const &eval_ddshape(xi_t const &)
	{
		return m_ddshape;
	}

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return &(domain_t::get_center());
	}

protected:
	/** \brief static shape function vector */
	static const shape_t m_shape;
	/** \brief static shape function gradient */
	static const dshape_t m_dshape;
	/** \brief static shape function second derivatives */
	static const ddshape_t m_ddshape;
};


/** constant line shape set */
typedef constant_shape_set<line_domain> line_0_shape_set;
/** constant triangle shape set */
typedef constant_shape_set<tria_domain> tria_0_shape_set;
/** constant quad shape set */
typedef constant_shape_set<quad_domain> quad_0_shape_set;
/** constant brick shape set */
typedef constant_shape_set<brick_domain> brick_0_shape_set;


/** \brief shape functions of the constant line shape set */
template <>
typename line_0_shape_set::shape_t
	const line_0_shape_set::m_shape = line_0_shape_set::shape_t::Ones();

/** \brief shape function derivatives of the constant line shape set */
template <>
typename line_0_shape_set::dshape_t
	const line_0_shape_set::m_dshape = line_0_shape_set::dshape_t::Zero();

/** \brief shape function second derivatives of the constant line shape set */
template <>
typename line_0_shape_set::ddshape_t
	const line_0_shape_set::m_ddshape = line_0_shape_set::ddshape_t::Zero();

/** \brief shape functions of the constant quad shape set */
template <>
typename tria_0_shape_set::shape_t
	const tria_0_shape_set::m_shape = tria_0_shape_set::shape_t::Ones();

/** \brief shape function derivatives of the constant quad shape set */
template <>
typename tria_0_shape_set::dshape_t
	const tria_0_shape_set::m_dshape = tria_0_shape_set::dshape_t::Zero();

/** \brief second shape function derivatives of the constant quad shape set */
template <>
typename tria_0_shape_set::ddshape_t
	const tria_0_shape_set::m_ddshape = tria_0_shape_set::ddshape_t::Zero();

/** \brief shape functions of the constant tria shape set */
template <>
typename quad_0_shape_set::shape_t
	const quad_0_shape_set::m_shape = quad_0_shape_set::shape_t::Ones();

/** \brief shape function derivatives of the constant tria shape set */
template <>
typename quad_0_shape_set::dshape_t
	const quad_0_shape_set::m_dshape = quad_0_shape_set::dshape_t::Zero();

/** \brief second shape function derivatives of the constant tria shape set */
template <>
typename quad_0_shape_set::ddshape_t
	const quad_0_shape_set::m_ddshape = quad_0_shape_set::ddshape_t::Zero();

/** \brief shape functions of the constant brick shape set */
template <>
typename brick_0_shape_set::shape_t
	const brick_0_shape_set::m_shape = brick_0_shape_set::shape_t::Ones();

/** \brief shape function derivatives of the constant brick shape set */
template <>
typename brick_0_shape_set::dshape_t
	const brick_0_shape_set::m_dshape = brick_0_shape_set::dshape_t::Zero();

/** \brief second shape function derivatives of the constant brick shape set */
template <>
typename brick_0_shape_set::ddshape_t
	const brick_0_shape_set::m_ddshape = brick_0_shape_set::ddshape_t::Zero();



/** \brief Common traits for all isoparametric shape sets
* \tparam Domain the domain class of the shape set
*/
template <class Domain>
struct isoparam_shape_set_traits_base
{
	/** \brief the domain type */
	typedef Domain domain_t;

	enum {
		num_nodes = domain_t::num_corners,
		polynomial_order = 1
	};
};


// Forward declaration
template <class Domain>
class isoparam_shape_set;

/** \brief Traits of the isoparametric line shape set */
template <>
struct shape_set_traits<isoparam_shape_set<line_domain> >
	: isoparam_shape_set_traits_base<line_domain>
{
	/** \brief highest power of local variables in the jacobian */
	enum { jacobian_order = 0 };
};

/** \brief Traits of the isoparametric triangle shape set */
template <>
struct shape_set_traits<isoparam_shape_set<tria_domain> >
	: isoparam_shape_set_traits_base<tria_domain>
{
	/** \brief highest power of local variables in the jacobian */
	enum { jacobian_order = 0 };
};

/** \brief Traits of the isoparametric quadrangle shape set */
template <>
struct shape_set_traits<isoparam_shape_set<quad_domain> >
	: isoparam_shape_set_traits_base<quad_domain>
{
	/** \brief highest power of local variables in the jacobian */
	enum { jacobian_order = 1 };
};

/** \brief Traits of the isoparametric brick shape set */
template <>
struct shape_set_traits<isoparam_shape_set<brick_domain> >
	: isoparam_shape_set_traits_base<brick_domain>
{
	/** \brief highest power of local variables in the jacobian */
	enum { jacobian_order = 1 };
};

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
	/** \brief type of a shape vector */
	typedef typename base_t::shape_t shape_t;
	/** \brief type of a shape gradient matrix */
	typedef typename base_t::dshape_t dshape_t;
	/** \brief type of a second shape derivative matrix */
	typedef typename base_t::ddshape_t ddshape_t;

	/**
	* \brief return shape functions
	* \param [in] xi the location in the base domain
	* \return the shape functions
	*/
	static shape_t eval_shape(xi_t const &xi);

	/**
	* \brief return shape functions derivatives
	* \param [in] xi the location in the base domain
	* \return the shape function derivatives
	*/
	static dshape_t eval_dshape(xi_t const &xi);

	/**
	* \brief return second shape function derivatives
	* \param [in] xi the location in the base domain
	* \return second shape function derivatives
	*/
	static ddshape_t eval_ddshape(xi_t const &xi);

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

/**
* \brief linear 2-noded line shape functions
* \param [in] xi domain location vector
* \return shape function vector set
* \details The shape functions are
*
* \f$L_1(\xi) = (1-\xi)/2 \\ L_2(\xi) = (1+\xi)/2 \f$
*/
template<>
inline typename line_1_shape_set::shape_t
	line_1_shape_set::eval_shape(typename line_1_shape_set::xi_t const &xi)
{
	shape_t L;
	L <<
		(1.0-xi[0])/2.0,
		(1.0+xi[0])/2.0;
	return L;
}

/**
* \brief linear 2-noded line shape function derivative matrix
* \return shape function derivative matrix
* \details The shape functions are
*
* \f$L'_1(\xi) = -1/2 \\ L'_2(\xi) = 1/2 \f$
*/
template<>
inline typename isoparam_shape_set<line_domain>::dshape_t
	line_1_shape_set::eval_dshape(typename line_1_shape_set::xi_t const &)
{
	dshape_t dL;
	dL <<
		-0.5,
		+0.5;
	return dL;
}


/**
* \brief linear 2-noded line shape function second derivative matrix
* \return shape function second derivative matrix
* \details The shape function derivatives are
*
* \f$L''_1(\xi) = 0 \f$
*/
template<>
inline typename isoparam_shape_set<line_domain>::ddshape_t
	line_1_shape_set::eval_ddshape(typename line_1_shape_set::xi_t const &)
{
	return ddshape_t::Zero();
}


/** linear tria shape set */
typedef isoparam_shape_set<tria_domain> tria_1_shape_set;

/**
* \brief linear 3-noded triangle shape functions
* \param [in] xi the domain variable vector
* \return the shape function vector
* \details The shape functions are
*
* \f$L_1(\xi, \eta) = 1-\xi-\eta \\ L_2(\xi, \eta) = \xi \\ L_3(\xi, \eta) = \eta \f$
*/
template<>
inline typename tria_1_shape_set::shape_t
	tria_1_shape_set::eval_shape(typename tria_1_shape_set::xi_t const &xi)
{
	shape_t L;
	L <<
		1.0-xi[0]-xi[1],
		xi[0],
		xi[1];
	return L;
}

/**
* \brief linear 3-noded tria elem shape function derivative matrix
* \return sape function derivative matrix
*/
template<>
inline typename tria_1_shape_set::dshape_t
	tria_1_shape_set::eval_dshape(typename tria_1_shape_set::xi_t const &)
{
	dshape_t dL;
	dL <<
		-1.0, -1.0,
		+1.0,  0.0,
		0.0, +1.0;
	return dL;
}

/**
* \brief linear 3-noded tria elem shape function second derivative matrix
* \return sape function second derivative matrix
*/
template<>
inline typename tria_1_shape_set::ddshape_t
	tria_1_shape_set::eval_ddshape(typename tria_1_shape_set::xi_t const &)
{
	return ddshape_t::Zero();
}



/** linear quad shape set */
typedef isoparam_shape_set<quad_domain> quad_1_shape_set;

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
inline typename quad_1_shape_set::shape_t
	quad_1_shape_set::eval_shape(typename quad_1_shape_set::xi_t const &xi)
{
	shape_t L;
	L <<
		(1.0-xi[0])*(1.0-xi[1])/4.0,
		(1.0+xi[0])*(1.0-xi[1])/4.0,
		(1.0+xi[0])*(1.0+xi[1])/4.0,
		(1.0-xi[0])*(1.0+xi[1])/4.0;
	return L;
}

/**
* \brief linear 4-noded general quadrilater shape function derivative matrix
* \return shape function gradient matrix
*/
template<>
inline typename quad_1_shape_set::dshape_t
	quad_1_shape_set::eval_dshape(typename quad_1_shape_set::xi_t const &xi)
{
	dshape_t dL;
	dL <<
		-.25 * (1-xi[1]), (1.0-xi[0]) * -.25,
		+.25 * (1-xi[1]), (1.0+xi[0]) * -.25,
		+.25 * (1+xi[1]), (1.0+xi[0]) * +.25,
		-.25 * (1+xi[1]), (1.0-xi[0]) * +.25;
	return dL;
}

/**
* \brief linear 4-noded general quadrilater shape function second derivative matrix
* \return shape function second derivative matrix
*/
template<>
inline typename quad_1_shape_set::ddshape_t
	quad_1_shape_set::eval_ddshape(typename quad_1_shape_set::xi_t const &xi)
{
	ddshape_t ddL;
	ddL <<
		0.0, +.25, 0.0,
		0.0, -.25, 0.0,
		0.0, +.25, 0.0,
		0.0, -.25, 0.0;
	return ddL;
}



/** linear brick shape set */
typedef isoparam_shape_set<brick_domain> brick_1_shape_set;

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
inline typename brick_1_shape_set::shape_t
	brick_1_shape_set::eval_shape(typename brick_1_shape_set::xi_t const &xi)
{
	shape_t L;
	L <<
		(1.0-xi[0])*(1.0-xi[1])*(1.0-xi[2])/8.0,
		(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[2])/8.0,
		(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[2])/8.0,
		(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[2])/8.0,
		(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[2])/8.0,
		(1.0+xi[0])*(1.0-xi[1])*(1.0+xi[2])/8.0,
		(1.0+xi[0])*(1.0+xi[1])*(1.0+xi[2])/8.0,
		(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[2])/8.0;
	return L;
}

/**
* \brief linear 8-noded general brick shape function derivative matrix
* \param [in] xi the domain variable vector
* \return shape function gradient matrix
*/
template<>
inline typename brick_1_shape_set::dshape_t
	brick_1_shape_set::eval_dshape(typename brick_1_shape_set::xi_t const &xi)
{
	dshape_t dL;
	dL <<
		(-1.0)*(1.0-xi[1])*(1.0-xi[2])/8.0, (1.0-xi[0])*(-1.0)*(1.0-xi[2])/8.0, (1.0-xi[0])*(1.0-xi[1])*(-1.0)/8.0,
		(+1.0)*(1.0-xi[1])*(1.0-xi[2])/8.0, (1.0+xi[0])*(-1.0)*(1.0-xi[2])/8.0, (1.0+xi[0])*(1.0-xi[1])*(-1.0)/8.0,
		(+1.0)*(1.0+xi[1])*(1.0-xi[2])/8.0, (1.0+xi[0])*(+1.0)*(1.0-xi[2])/8.0, (1.0+xi[0])*(1.0+xi[1])*(-1.0)/8.0,
		(-1.0)*(1.0+xi[1])*(1.0-xi[2])/8.0, (1.0-xi[0])*(+1.0)*(1.0-xi[2])/8.0, (1.0-xi[0])*(1.0+xi[1])*(-1.0)/8.0,
		(-1.0)*(1.0-xi[1])*(1.0+xi[2])/8.0, (1.0-xi[0])*(-1.0)*(1.0+xi[2])/8.0, (1.0-xi[0])*(1.0-xi[1])*(+1.0)/8.0,
		(+1.0)*(1.0-xi[1])*(1.0+xi[2])/8.0, (1.0+xi[0])*(-1.0)*(1.0+xi[2])/8.0, (1.0+xi[0])*(1.0-xi[1])*(+1.0)/8.0,
		(+1.0)*(1.0+xi[1])*(1.0+xi[2])/8.0, (1.0+xi[0])*(+1.0)*(1.0+xi[2])/8.0, (1.0+xi[0])*(1.0+xi[1])*(+1.0)/8.0,
		(-1.0)*(1.0+xi[1])*(1.0+xi[2])/8.0, (1.0-xi[0])*(+1.0)*(1.0+xi[2])/8.0, (1.0-xi[0])*(1.0+xi[1])*(+1.0)/8.0;
	return dL;
}

/**
* \brief linear 8-noded general brick shape function second derivative matrix
* \param [in] xi the domain variable vector
* \return shape function second derivative matrix
*/
template<>
inline typename brick_1_shape_set::ddshape_t
	brick_1_shape_set::eval_ddshape(typename brick_1_shape_set::xi_t const &xi)
{
	ddshape_t dL;
	dL <<
		0,  1/8.0 - xi[2]/8.0,  1/8.0 - xi[1]/8.0, 0,  1/8.0 - xi[0]/8.0, 0,
		0,  xi[2]/8.0 - 1/8.0,  xi[1]/8.0 - 1/8.0, 0,  xi[0]/8.0 + 1/8.0, 0,
		0,  1/8.0 - xi[2]/8.0, -xi[1]/8.0 - 1/8.0, 0, -xi[0]/8.0 - 1/8.0, 0,
		0,  xi[2]/8.0 - 1/8.0,  xi[1]/8.0 + 1/8.0, 0,  xi[0]/8.0 - 1/8.0, 0,
		0,  xi[2]/8.0 + 1/8.0,  xi[1]/8.0 - 1/8.0, 0,  xi[0]/8.0 - 1/8.0, 0,
		0, -xi[2]/8.0 - 1/8.0,  1/8.0 - xi[1]/8.0, 0, -xi[0]/8.0 - 1/8.0, 0,
		0,  xi[2]/8.0 + 1/8.0,  xi[1]/8.0 + 1/8.0, 0,  xi[0]/8.0 + 1/8.0, 0,
		0, -xi[2]/8.0 - 1/8.0, -xi[1]/8.0 - 1/8.0, 0,  1/8.0 - xi[0]/8.0, 0;
	return dL;
}


// Forward declaration
class parallelogram_shape_set;

/** \brief Traits for parallelogram shapesets */
template<>
struct shape_set_traits<parallelogram_shape_set>
{
	/** \brief the domain type */
	typedef quad_domain domain_t;
	/** \brief number of nodes */
	enum {
		num_nodes = 3,
		polynomial_order = 1,
		jacobian_order = 0
	};
};

/**
* \brief linear 3-noded parallelogram shape function set
*/
class parallelogram_shape_set : public shape_set_base<parallelogram_shape_set>
{
public:
	/**
	* \brief linear 3-noded parallelogram shape functions
	* \param [in] xi the domain variable
	* \return the shape function vector
	* \details The shape functions are
	*
	* \f$L_1(\xi, \eta) = (-\xi-\eta)/2 \\ L_2(\xi, \eta) = (1+\xi)/2 \\ L_3(\xi, \eta) = (1+\eta)/2 \f$
	*/
	static shape_t eval_shape(xi_t const &xi)
	{
		shape_t L;
		L <<
			(-xi[0]-xi[1])/2.0,
			(1.0+xi[0])/2.0,
			(1.0+xi[1])/2.0;
		return L;
	}

	/**
	* \brief linear 3-noded parallelogram shape function derivatives
	* \return the shape function gradient matrix
	*/
	static dshape_t eval_dshape(xi_t const &)
	{
		dshape_t dL;
		dL <<
			-.5, -.5,
			+.5,  .0,
			.0, +.5;
		return dL;
	}

	/**
	* \brief linear 3-noded parallelogram second shape function derivatives
	* \return the shape function second derivative matrix
	*/
	static ddshape_t eval_ddshape(xi_t const &)
	{
		return ddshape_t::Zero();
	}

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return domain_t::get_corners();
	}
};



// Forward declaration
class line_2_shape_set;

/** \brief Traits for quadratic line shapesets */
template<>
struct shape_set_traits<line_2_shape_set>
{
	/** \brief the domain type */
	typedef line_domain domain_t;
	/** \brief number of nodes */
	enum {
		num_nodes = 3,
		polynomial_order = 2,
		jacobian_order = 1
	};
};

/**
* \brief quadratic 3-noded line shape function set
*/
class line_2_shape_set : public shape_set_base<line_2_shape_set>
{
public:
	/**
	* \brief quadratic 3-noded line shape functions
	* \param [in] _xi the domain variable
	* \return the shape function vector
	*/
	static shape_t eval_shape(xi_t const &_xi)
	{
		scalar_t xi = _xi[0];
		shape_t L;
		L <<
		 -xi*(1.0-xi)/2.0,
		 1.0-xi*xi,
		 xi*(1.0+xi)/2.0;
		return L;
	}

	/**
	* \brief quadratic 3-noded line shape function derivatives
	* \param [in] _xi the domain variable
	* \return the shape function gradient matrix
	*/
	static dshape_t eval_dshape(xi_t const & _xi)
	{
		scalar_t xi = _xi[0];
		dshape_t dL;
		dL <<
			xi - 0.5,
			-2.0*xi,
			xi + 0.5;
		return dL;
	}

	/**
	* \brief quadratic 3-noded line shape function second derivatives
	* \return the shape function second derivative matrix
	*/
	static ddshape_t eval_ddshape(xi_t const &)
	{
		return ddshape_t(1.0, -2.0, 1.0);
	}

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

int const line_2_shape_set::m_domain_indices[line_2_shape_set::num_nodes] = {0, -1, 1};



// Forward declaration
class tria_2_shape_set;

/** \brief Traits for quadratic tria shapesets */
template<>
struct shape_set_traits<tria_2_shape_set>
{
	/** \brief the domain type */
	typedef tria_domain domain_t;
	/** \brief number of nodes */
	enum {
		num_nodes = 6,
		polynomial_order = 2,
		jacobian_order = 2
	};
};

/**
* \brief quadratic 6-noded tria shape function set
*/
class tria_2_shape_set : public shape_set_base<tria_2_shape_set>
{
public:
	/**
	* \brief quadratic 6-noded tria shape functions
	* \param [in] _xi the domain variable
	* \return the shape function vector
	*/
	static shape_t eval_shape(xi_t const &_xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1];
		shape_t L;
		L <<
			(eta+xi-1)*(2*eta+2*xi-1),
			-4*xi*(eta+xi-1),
			xi*(2*xi-1),
			4*eta*xi,
			eta*(2*eta-1),
			-4*eta*(eta+xi-1);
		return L;
	}

	/**
	* \brief quadratic 6-noded tria shape function derivatives
	* \param [in] _xi the domain variable
	* \return the shape function gradient matrix
	*/
	static dshape_t eval_dshape(xi_t const & _xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1];
		dshape_t dL;
		dL <<
			4*eta+4*xi-3, 4*eta+4*xi-3,
			4-8*xi-4*eta, -4*xi,
			4*xi-1,       0,
			4*eta,        4*xi,
			0,            4*eta-1,
			-4*eta,       4-4*xi-8*eta;
		return dL;
	}

	/**
	 * \brief quadratic 6-noded tria shape function second derivatives
	 * \return the shape function second derivative matrix
	 */
	static ddshape_t eval_ddshape(xi_t const &)
	{
		ddshape_t ddL;
		ddL <<
			 4,  4,  4,
			-8, -4,  0,
			 4,  0,  0,
			 0,  4,  0,
			 0,  0,  4,
			 0, -4, -8;
		return ddL;
	}

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
	{0, -1, 1, -1, 2, -1};



// Forward declaration
class quad_2_shape_set;

/** \brief Traits for quadratic quad shapesets */
template<>
struct shape_set_traits<quad_2_shape_set>
{
	typedef quad_domain domain_t;	/**< \brief the domain type */
	enum {
		num_nodes = 9,
		polynomial_order = 2,
		jacobian_order = 3
	};
};

/**
* \brief quadratic 9-noded quad shape function set
*/
class quad_2_shape_set : public shape_set_base<quad_2_shape_set>
{
public:
	/**
	* \brief quadratic 9-noded quad shape functions
	* \param [in] _xi the domain variable
	* \return the shape function vector
	*/
	static shape_t eval_shape(xi_t const &_xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1];
		scalar_t _1mxi = 1-xi, _1pxi = 1+xi;
		scalar_t _1meta = 1-eta, _1peta = 1+eta;
		shape_t L;
		L <<
			_1mxi*xi * _1meta*eta/4.0,
			_1mxi*_1pxi * _1meta*(-eta)/2.0,
			_1pxi*xi * _1meta*(-eta)/4.0,
			_1pxi*xi * _1meta*_1peta/2.0,
			_1pxi*xi * _1peta*eta/4.0,
			_1mxi*_1pxi * _1peta*eta/2.0,
			_1mxi*(-xi) * _1peta*eta/4.0,
			_1mxi*(-xi) * _1meta*_1peta/2.0,
			_1mxi*_1pxi * _1meta*_1peta;
		return L;
	}

	/**
	 * \brief quadratic 9-noded quad shape function derivatives
	 * \param [in] _xi the domain variable
	 * \return the shape function gradient matrix
	 */
	static dshape_t eval_dshape(xi_t const & _xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1], xi2 = xi*xi, eta2 = eta*eta;
		dshape_t dL;
		dL <<
            eta*(2.0*xi-1.0)*(eta-1.0)/4.0, xi*(2.0*eta-1.0)*(xi-1.0)/4.0,
            -xi*eta*(eta-1.0),              -(xi2-1.0)*(2.0*eta-1.0)/2.0,
            eta*(2.0*xi+1.0)*(eta-1.0)/4.0, xi*(2.0*eta-1.0)*(xi+1.0)/4.0,
            -(2.0*xi+1.0)*(eta2-1.0)/2.0, -xi*eta*(xi+1.0),
            eta*(2.0*xi+1.0)*(eta+1.0)/4.0, xi*(2.0*eta+1.0)*(xi+1.0)/4.0,
            -xi*eta*(eta+1.0),              -(xi2-1.0)*(2.0*eta+1.0)/2.0,
            eta*(2.0*xi-1.0)*(eta+1.0)/4.0, xi*(2.0*eta+1.0)*(xi-1.0)/4.0,
            -(2.0*xi-1.0)*(eta2-1.0)/2.0,   -xi*eta*(xi-1.0),
            2.0*xi*(eta2-1.0),              2.0*eta*(xi2-1.0);
		return dL;
	}

	/**
	 * \brief quadratic 9-noded quad shape function second derivatives
	 * \param [in] _xi the domain variable
	 * \return the shape function second derivative matrix
	 */
	static ddshape_t eval_ddshape(xi_t const & _xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1], xi2 = xi*xi, eta2 = eta*eta;
		ddshape_t ddL;
		ddL <<
			eta2/2.0 - eta/2.0, eta*xi - xi/2.0 - eta/2.0 + 1.0/4.0, xi2/2.0 - xi/2.0,
			- eta2 + eta,   xi - 2*eta*xi,               1.0 - xi2,
			eta2/2.0 - eta/2.0, eta/2.0 - xi/2.0 + eta*xi - 1.0/4.0, xi2/2.0 + xi/2.0,
			1.0 - eta2,       - eta - 2*eta*xi,            - xi2 - xi,
			eta2/2.0 + eta/2.0, eta/2.0 + xi/2.0 + eta*xi + 1.0/4.0, xi2/2.0 + xi/2.0,
			- eta2 - eta,   - xi - 2*eta*xi,             1.0 - xi2,
			eta2/2.0 + eta/2.0, xi/2.0 - eta/2.0 + eta*xi - 1.0/4.0, xi2/2.0 - xi/2.0,
			1.0 - eta2,       eta - 2*eta*xi,              - xi2 + xi,
			2*eta2 - 2.0,     4*eta*xi,                    2*xi2 - 2.0;
   		return ddL;
	}

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
		quad_2_shape_set::xi_t(-1.0,-1.0),
		quad_2_shape_set::xi_t( 0.0,-1.0),
		quad_2_shape_set::xi_t(+1.0,-1.0),
		quad_2_shape_set::xi_t(+1.0, 0.0),
		quad_2_shape_set::xi_t(+1.0,+1.0),
		quad_2_shape_set::xi_t( 0.0,+1.0),
		quad_2_shape_set::xi_t(-1.0,+1.0),
		quad_2_shape_set::xi_t(-1.0, 0.0),
		quad_2_shape_set::xi_t( 0.0, 0.0)
};

int const quad_2_shape_set::m_domain_corners[quad_2_shape_set::num_nodes] =
	{0, -1, 1, -1, 2, -1, 3, -1, -1};



#endif // SHAPESET_HPP_INCLUDED

