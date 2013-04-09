/**
 * \file shapeset.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief Definition of various shape function sets
 */
#ifndef SHAPESET_HPP_INCLUDED
#define SHAPESET_HPP_INCLUDED

#include "domain.hpp"
#include "../tmp/bool.hpp"

/**
 * \brief Traits for shapesets
 */
template <class Dervied>
struct shape_set_traits;

/**
 * \brief Shapeset base class for CRTP
 * \tparam Derived The derived shapeset class
 */
template <class Derived>
class shape_set_base
{
public:
	typedef typename shape_set_traits<Derived>::domain_t domain_t;			/** \brief Domain */
	static unsigned const num_nodes = shape_set_traits<Derived>::num_nodes;	/** \brief Number of nodes */

	/** \brief type of the local coordinate */
	typedef typename domain_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector */
	typedef Eigen::Matrix<double, num_nodes, 1> shape_t;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix */
	typedef Eigen::Matrix<double, num_nodes, domain_t::dimension> dshape_t;

public:
	/**
	 *\brief shape function vector \fL_i(\xi)\f$
	 */
	static shape_t eval_shape(xi_t const &xi)
	{
		return Derived::eval_shape(xi);
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 */
	static dshape_t eval_dshape(xi_t const &xi)
	{
		return Derived::eval_dshape(xi);
	}

	/**
	 * \brief begin iterator of corner nodes
	 */
	static xi_t const* corner_begin(void)
	{
		return Derived::corner_begin_impl();
	}

	/**
	 * \brief end iterator of corner nodes
	 */
	static xi_t const* corner_end(void)
	{
		return Derived::corner_end_impl();
	}
};

// Forward declaration
template <class Domain>
class constant_shape_set;

/**
 * \brief Traits for constant shape sets
 */
template <class Domain>
struct shape_set_traits<constant_shape_set<Domain> >
{
	typedef Domain domain_t;
	static unsigned const num_nodes = 1;
};

/**
 * \brief Constant shape sets
 */
template <class Domain>
class constant_shape_set : public shape_set_base<constant_shape_set<Domain> >
{
public:
	typedef shape_set_base<constant_shape_set<Domain> > base_t;
	typedef typename base_t::domain_t domain_t;
	typedef typename base_t::xi_t xi_t;
	typedef typename base_t::shape_t shape_t;
	typedef typename base_t::dshape_t dshape_t;

public:
	/**
	 * \brief Constant shape functions
	 * \details The shape functions are
	 *
	 * \f$L_1(\xi) = 1 \f$
	 */
	static shape_t eval_shape(xi_t const &)
	{
		return shape_t::Ones();
	}

	/**
	 * \brief Derivatives of constant shape functions
	 * \details The derivatives are
	 *
	 * \f$\nabla L_1(\xi) = 0\f$
	 */
	static dshape_t eval_dshape(xi_t const &xi)
	{
		return shape_t::Zero();
	}

	static xi_t const *corner_begin_impl(void)
	{
		return &(domain_t::get_center());
	}

	static xi_t const *corner_end_impl(void)
	{
		return &(domain_t::get_center()) + 1;
	}
};

// Forward declaration
template <class Domain>
class isoparam_shape_set;

/**
 * \brief Traits for isoparametric shape sets
 */
template <class Domain>
struct shape_set_traits<isoparam_shape_set<Domain> >
{
	typedef Domain domain_t;
	static unsigned const num_nodes = domain_t::id; //Note: ID has actually nothing to do with id
};

/**
 * \brief Isoparametric shape sets
 */
template <class Domain>
class isoparam_shape_set : public shape_set_base<isoparam_shape_set<Domain> >
{
public:
	typedef Domain domain_t;

	typedef shape_set_base<isoparam_shape_set<Domain> > base_t;
	typedef typename base_t::xi_t xi_t;
	typedef typename base_t::shape_t shape_t;
	typedef typename base_t::dshape_t dshape_t;

	static shape_t eval_shape(xi_t const &xi);
	static dshape_t eval_dshape(xi_t const &xi);

	static xi_t const * corner_begin_impl(void)
	{
		return domain_t::get_corners();
	}

	static xi_t const * corner_end_impl(void)
	{
		return domain_t::get_corners() + domain_t::id;	//Note: ID has actually nothing to do with id, its the number of corners
	}
};

/**
 * \brief linear 2-noded line shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (1-\xi)/2 \\ L_2(\xi) = (1+\xi)/2 \f$
 */
template<>
typename isoparam_shape_set<line_domain>::shape_t
isoparam_shape_set<line_domain>::eval_shape(typename isoparam_shape_set<line_domain>::xi_t const &xi)
{
	shape_t L;
	L <<
		(1.0-xi[0])/2.0,
		(1.0+xi[0])/2.0;
	return L;
}

/**
 * \brief linear 2-noded line shape function derivative matrix
 */
template<>
typename isoparam_shape_set<line_domain>::dshape_t
isoparam_shape_set<line_domain>::eval_dshape(typename isoparam_shape_set<line_domain>::xi_t const &)
{
	dshape_t dL;
	dL <<
		-0.5,
		+0.5;
	return dL;
}


/**
 * \brief linear 3-noded triangle shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = 1-\xi-\eta \\ L_2(\xi) = \xi \\ L_3(\xi) = \eta \f$
 */
template<>
typename isoparam_shape_set<tria_domain>::shape_t
isoparam_shape_set<tria_domain>::eval_shape(typename isoparam_shape_set<tria_domain>::xi_t const &xi)
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
 */
template<>
typename isoparam_shape_set<tria_domain>::dshape_t
isoparam_shape_set<tria_domain>::eval_dshape(typename isoparam_shape_set<tria_domain>::xi_t const &)
{
	dshape_t dL;
	dL <<
		-1.0, -1.0,
		+1.0,  0.0,
		 0.0, +1.0;
	return dL;
}

/**
 * \brief linear 4-noded general quadrilateral shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (1-\xi)(1-\eta)/4 \\ L_2(\xi) = (1+\xi)(1-\eta)/4 \\ L_3(\xi) = (1+\xi)(1+\eta)/4 \\ L_4(\xi) = (1-\xi)(1+\eta)/4 \f$
 */
template<>
typename isoparam_shape_set<quad_domain>::shape_t
isoparam_shape_set<quad_domain>::eval_shape(typename isoparam_shape_set<quad_domain>::xi_t const &xi)
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
 */
template<>
typename isoparam_shape_set<quad_domain>::dshape_t
isoparam_shape_set<quad_domain>::eval_dshape(typename isoparam_shape_set<quad_domain>::xi_t const &xi)
{
	dshape_t dL;
	dL <<
		-.25 * (1-xi[1]), (1.0-xi[0]) * -.25,
		+.25 * (1-xi[1]), (1.0+xi[0]) * -.25,
		+.25 * (1+xi[1]), (1.0+xi[0]) * +.25,
		-.25 * (1+xi[1]), (1.0-xi[0]) * +.25;
	return dL;
}

// Forward declaration
class parallelogram_shape_set;

template<>
struct shape_set_traits<parallelogram_shape_set>
{
	typedef quad_domain domain_t;
	static unsigned const num_nodes = 3; //Note: ID has actually nothing to do with id
};

/**
 * \brief linear 3-noded parallelogram shape function set
 */
class parallelogram_shape_set : public shape_set_base<parallelogram_shape_set>
{
public:
	typedef shape_set_base<parallelogram_shape_set> base_t;
	typedef base_t::domain_t domain_t;
	typedef base_t::xi_t xi_t;
	typedef base_t::shape_t shape_t;
	typedef base_t::dshape_t dshape_t;
public:
	/**
	 * \brief linear 3-noded parallelogram shape functions
	 * \details The shape functions are
	 *
	 * \f$L_1(\xi) = (-\xi-\eta)/2 \\ L_2(\xi) = (1+\xi)/2 \\ L_3(\xi) = (1+\eta)/2 \f$
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
	static dshape_t eval_dshape(xi_t const &)
	{
		dshape_t dL;
		dL <<
			-.5, -.5,
			+.5,  .0,
			.0, +.5;
		return dL;
	}

	static xi_t const * corner_begin_impl(void)
	{
		return domain_t::get_corners();
	}

	static xi_t const * corner_end_impl(void)
	{
		return domain_t::get_corners() + domain_t::id;	//Note: ID has actually nothing to do with id, its the number of corners
	}
};

/**
 * \brief Type definitions for further usage
 */
typedef constant_shape_set<line_domain> line_0_shape_set;
typedef constant_shape_set<tria_domain> tria_0_shape_set;
typedef constant_shape_set<quad_domain> quad_0_shape_set;

typedef isoparam_shape_set<line_domain> line_1_shape_set;
typedef isoparam_shape_set<tria_domain> tria_1_shape_set;
typedef isoparam_shape_set<quad_domain> quad_1_shape_set;


#endif // SHAPESET_HPP_INCLUDED

/*
// TODO - use x_t instead of nDim template parameter
template <class shape_set_from, class shape_set_to, unsigned nDim>
struct shape_set_converter;

template <class Shape, unsigned nDim>
struct shape_set_converter<Shape, Shape, nDim>
{
	typedef Shape from_set;
	typedef Shape to_set;

	typedef Eigen::Matrix<double, from_set::num_nodes, nDim> from_coords_t;
	typedef Eigen::Matrix<double, to_set::num_nodes, nDim> to_coords_t;

public:
	static bool eval(from_coords_t const &coords)
	{
		to_coords = coords;
		return true;
	}

	static to_coords_t const &get_coords(void)
	{
		return to_coords;
	}

protected:
	static to_coords_t to_coords;
};

template <class Shape, unsigned nDim>
typename shape_set_converter<Shape, Shape, nDim>::to_coords_t
	shape_set_converter<Shape, Shape, nDim>::to_coords;



template <unsigned nDim>
struct shape_set_converter<quad_1_shape_set, parallelogram_shape_set, nDim>
{
	typedef quad_1_shape_set from_set;
	typedef parallelogram_shape_set to_set;

	typedef Eigen::Matrix<double, from_set::num_nodes, nDim> from_coords_t;
	typedef Eigen::Matrix<double, to_set::num_nodes, nDim> to_coords_t;

public:
	static bool eval(from_coords_t const &coords)
	{

		to_coords.row(0) = coords.row(0);
		to_coords.row(1) = coords.row(1);
		to_coords.row(2) = coords.row(3);

		to_set::xi_t xi;
		xi << 1.0, 1.0;

		return (to_set::eval_L(xi).transpose() * to_coords - coords.row(2)).norm() < 1e-3;
	}

	static to_coords_t const &get_coords(void)
	{
		return to_coords;
	}

protected:
	static to_coords_t to_coords;
};

template <unsigned nDim>
typename shape_set_converter<quad_1_shape_set, parallelogram_shape_set, nDim>::to_coords_t
	shape_set_converter<quad_1_shape_set, parallelogram_shape_set, nDim>::to_coords;

*/

