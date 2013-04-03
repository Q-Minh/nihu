/**
 * \file shapeset.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief Definition of class shape_set and its specialisations
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
class shape_set
{
public:
	/** \brief the shape function's domain type */
	typedef Domain domain_t;
	/** \brief number of shape functions in the set */
	static unsigned const num_nodes = NumNodes;

	/** \brief type of the shape functions' independent variable \f$\xi\f$ */
	typedef typename domain_t::xi_t xi_t;
	/** \brief type of corner storage */
	typedef xi_t corners_t[NumNodes];
	/** \brief type of an \f$L(\xi)\f$ vector */
	typedef Eigen::Matrix<double, num_nodes, 1> shape_t;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix */
	typedef Eigen::Matrix<double, num_nodes, domain_t::dimension> dshape_t;
	
	static xi_t const *corner_begin(void)
	{
		return m_corners;
	}

	static xi_t const *corner_end(void)
	{
		return m_corners + num_nodes;
	}
protected:
	static corners_t const m_corners;
};


/**
 * \brief linear 2-noded line shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (1-\xi)/2 \\ L_2(\xi) = (1+\xi)/2 \f$
 */
class line_1_shape_set : public shape_set<line_domain, 2>
{
public:
	/**
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static shape_t eval_shape(xi_t const &xi)
	{
		shape_t L;
		L <<
			(1.0-xi[0])/2.0,
			(1.0+xi[0])/2.0;
		return L;
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 */
	static dshape_t eval_dshape(xi_t const &)
	{
		dshape_t dL;
		dL <<
			-0.5,
			+0.5;
		return dL;
	}
};

template <>
shape_set<line_domain,2>::corners_t const shape_set<line_domain,2>::m_corners = {line_1_shape_set::xi_t(-1.0), line_1_shape_set::xi_t(1.0)};


/**
 * \brief linear 3-noded triangle shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = 1-\xi-\eta \\ L_2(\xi) = \xi \\ L_3(\xi) = \eta \f$
 */
class tria_1_shape_set : public shape_set<tria_domain, 3>
{
public:
	/**
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static shape_t eval_shape(xi_t const &xi)
	{
		shape_t L;
		L <<
			1.0-xi[0]-xi[1],
			xi[0],
			xi[1];
		return L;
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \brief return shape function gradient matrix
	 */
	static dshape_t eval_dshape(xi_t const &)
	{
		dshape_t dL;
		dL <<
			-1.0, -1.0,
			+1.0,  0.0,
			 0.0, +1.0;
		return dL;
	}
};

template<>
shape_set<tria_domain, 3>::corners_t const shape_set<tria_domain, 3>::m_corners =  {tria_1_shape_set::xi_t(0.0,0.0), tria_1_shape_set::xi_t(1.0,0.0), tria_1_shape_set::xi_t(0.0,1.0)};


/**
 * \brief linear 3-noded parallelogram shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (-\xi-\eta)/2 \\ L_2(\xi) = (1+\xi)/2 \\ L_3(\xi) = (1+\eta)/2 \f$
 */
class parallelogram_shape_set : public shape_set<quad_domain, 3>
{
public:
	/**
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
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
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static dshape_t eval_dshape(xi_t const &xi)
	{
		dshape_t dL;
		dL <<
			-.5, -.5,
			+.5,  .0,
			 .0, +.5;
		return dL;
	}
};

template<>
shape_set<quad_domain, 3>::corners_t const shape_set<quad_domain, 3>::m_corners =  {parallelogram_shape_set::xi_t(-1.0,-1.0), parallelogram_shape_set::xi_t(1.0,-1.0), parallelogram_shape_set::xi_t(1.0,1.0)};


/**
 * \brief linear 4-noded general quadrilateral shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = (1-\xi)(1-\eta)/4 \\ L_2(\xi) = (1+\xi)(1-\eta)/4 \\ L_3(\xi) = (1+\xi)(1+\eta)/4 \\ L_4(\xi) = (1-\xi)(1+\eta)/4 \f$
 */
class quad_1_shape_set : public shape_set<quad_domain, 4>
{
public:
	/**
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static shape_t eval_shape(xi_t const &xi)
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
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static dshape_t eval_dshape(xi_t const &xi)
	{
		dshape_t dL;
		dL <<
			-.25 * (1-xi[1]), (1.0-xi[0]) * -.25,
			+.25 * (1-xi[1]), (1.0+xi[0]) * -.25,
			+.25 * (1+xi[1]), (1.0+xi[0]) * +.25,
			-.25 * (1+xi[1]), (1.0-xi[0]) * +.25;
		return dL;
	}
};

template<>
shape_set<quad_domain, 4>::corners_t const shape_set<quad_domain, 4>::m_corners =  {quad_1_shape_set::xi_t(-1.0,-1.0), quad_1_shape_set::xi_t(1.0,-1.0),
	quad_1_shape_set::xi_t(1.0,1.0), quad_1_shape_set::xi_t(-1.0,1.0)};


/**
 * \brief constant shape function set
 * \tparam Domain the domain over which the shape function set is defined
 */
template <class Domain>
class constant_shape_set : public shape_set<Domain, 1>
{
public:
	typedef Domain domain_t;	/**< \brief template argument as nested type */
	
	typedef shape_set<domain_t, 1> base;	/**< \brief the shape set's type */
	typedef typename base::shape_t shape_t;	/**< \brief the L type */
	typedef typename base::dshape_t dshape_t;	/**< \brief the dL type */
	typedef typename base::xi_t xi_t;	/**< \brief the xi type */

	typedef typename base::corners_t corners_t;
	
	/**
	 * \brief shape function vector \f$L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static shape_t eval_shape(xi_t const &xi)
	{
		return shape_set<domain_t, 1>::shape_t::Ones();
	}

	/**
	 * \brief shape function gradient matrix \f$\nabla L_i(\xi)\f$
	 * \param \xi independent variable \f$\xi\f$
	 */
	static dshape_t eval_dshape(xi_t const &xi)
	{
		return dshape_t::Zero();
	}
	/*
	static xi_t const *corner_begin(void)
	{
		return &(domain_t::get_center());
	}

	static xi_t const *corner_end(void)
	{
		return &(domain_t::get_center())+1;
	}
	*/
};


/**
 * \brief constant 3-noded triangular shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = 1 \f$
 */
typedef constant_shape_set<tria_domain> tria_0_shape_set;

template<>
shape_set<tria_domain, 1>::corners_t const shape_set<tria_domain, 1>::m_corners = {tria_domain::get_center()};

/**
 * \brief constant 4-noded quadrilateral shape functions
 * \details The shape functions are
 *
 * \f$L_1(\xi) = 1 \f$
 */
typedef constant_shape_set<quad_domain> quad_0_shape_set;

template<>
shape_set<quad_domain, 1>::corners_t const shape_set<quad_domain, 1>::m_corners = {quad_domain::get_center()};


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

#endif

