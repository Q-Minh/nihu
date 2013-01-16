/**
 * \file shapeset.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief Definition of class ShapeSet and its specialisations
 */
#ifndef SHAPESET_HPP_INCLUDED
#define SHAPESET_HPP_INCLUDED

#include "domain.hpp"
#include "../tmp/bool.hpp"

/**
 * \brief Shape function set
 * \tparam Domain the domain of the the shape function
 * \tparam NumNodes the number of shape functions in the set
 * \details This abstract class is used as the base class of all shape function sets.
 * The class defines some convenient types that are used in the derived definition.
 */
template <class Domain, unsigned NumNodes>
class ShapeSet
{
public:
	/** \brief the shape function's domain type */
	typedef Domain domain_t;
	/** \brief number of shape functions in the set */
	static unsigned const num_nodes = NumNodes;

	/** \brief type of the shape functions' independent variable \f$\xi\f$ */
	typedef typename domain_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector */
	typedef Matrix<double, num_nodes, 1> L_t;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix */
	typedef Matrix<double, num_nodes, domain_t::dimension> dL_t;
};


/**
 * \brief linear 2-noded line shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (1-\xi)/2 \\ L_2(\xi) = (1+\xi)/2 \f$
 */
class line_1_shape_set : public ShapeSet<line_domain, 2>
{
public:
	/** 
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static L_t eval_L(xi_t const &xi)
	{
		L_t L;
		L <<
			(1.0-xi[0])/2.0,
			(1.0+xi[0])/2.0;
		return L;
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static dL_t eval_dL(xi_t const &)
	{
		dL_t dL;
		dL <<
			-0.5,
			+0.5;
		return dL;
	}
};


/**
 * \brief linear 3-noded triangle shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = 1-\xi-\eta \\ L_2(\xi) = \xi \\ L_3(\xi) = \eta \f$
 */
class tria_1_shape_set : public ShapeSet<tria_domain, 3>
{
public:
	/** 
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static L_t eval_L(xi_t const &xi)
	{
		L_t L;
		L <<
			1.0-xi[0]-xi[1],
			xi[0],
			xi[1];
		return L;
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static dL_t eval_dL(xi_t const &)
	{
		dL_t dL;
		dL <<
			-1.0, -1.0,
			+1.0,  0.0,
			 0.0, +1.0;
		return dL;
	}
};


/**
 * \brief linear 3-noded parallelogram shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (-\xi-\eta)/2 \\ L_2(\xi) = (1+\xi)/2 \\ L_3(\xi) = (1+\eta)/2 \f$
 */
class parallelogram_shape_set : public ShapeSet<quad_domain, 3>
{
public:
	/** 
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static L_t eval_L(xi_t const &xi)
	{
		L_t L;
		L <<
			(-xi[0]-xi[1])/2.0,
			(1.0+xi[0])/2.0,
			(1.0+xi[1])/2.0;
		return L;
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static dL_t eval_dL(xi_t const &xi)
	{
		dL_t dL;
		dL <<
			-.5, -.5,
			+.5,  .0,
			 .0, +.5;
		return dL;
	}
};


/**
 * \brief linear 4-noded general quadrilateral shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (1-\xi)(1-\eta)/4 \\ L_2(\xi) = (1+\xi)(1-\eta)/4 \\ L_3(\xi) = (1+\xi)(1+\eta)/4 \\ L_4(\xi) = (1-\xi)(1+\eta)/4 \f$
 */
class quad_1_shape_set : public ShapeSet<quad_domain, 4>
{
public:
	/** 
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static L_t eval_L(xi_t const &xi)
	{
		L_t L;
		L <<
			(1.0-xi[0])*(1.0-xi[1])/4.0,
			(1.0+xi[0])*(1.0-xi[1])/4.0,
			(1.0+xi[0])*(1.0+xi[1])/4.0,
			(1.0-xi[0])*(1.0+xi[1])/4.0;
		return L;
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static dL_t eval_dL(xi_t const &xi)
	{
		dL_t dL;
		dL <<
			-.25 * (1-xi[1]), (1.0-xi[0]) * -.25,
			+.25 * (1-xi[1]), (1.0+xi[0]) * -.25,
			+.25 * (1+xi[1]), (1.0+xi[0]) * +.25,
			-.25 * (1+xi[1]), (1.0-xi[0]) * +.25;
		return dL;
	}
};

/**
 * \brief metafunction returns true if the argument is a linear shape function set
 * \tparam ShapeSet the shape function set class
 */
template <class ShapeSet>
struct is_linear : false_ {};
template <> struct is_linear<tria_1_shape_set> : true_ {};
template <> struct is_linear<parallelogram_shape_set> : true_ {};


template <class shape_set_from, class shape_set_to>
struct shape_set_converter;

template <>
struct shape_set_converter<quad_1_shape_set, parallelogram_shape_set>
{
	typedef quad_1_shape_set from_set;
	typedef parallelogram_shape_set to_set;
	
	template <int dim>
	static bool eval(Matrix<double, 4, dim> const &coords)
	{
		Matrix<double, 3, dim> c;
		c.row(0) = coords.row(0);
		c.row(1) = coords.row(1);
		c.row(2) = coords.row(3);
		
		to_set::xi_t xi;
		xi << 1.0, 1.0;
		
		return (
			to_set::eval_L(xi).transpose() * c - coords.row(2)
		).norm() < 1e-3;
	}
};

#endif

