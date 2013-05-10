/**
* \file shapeset.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief Definition of various shape function sets
*/
#ifndef SHAPESET_HPP_INCLUDED
#define SHAPESET_HPP_INCLUDED

#include "domain.hpp"
#include "../tmp/bool.hpp"

/** \brief Traits for shapesets */
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
	typedef shape_set_traits<Derived> traits_t;	/**< \brief the traits class type */
	typedef typename traits_t::domain_t domain_t;			/**< \brief Domain */
	static unsigned const num_nodes = traits_t::num_nodes;	/**< \brief Number of nodes */

	/** \brief scalar type inherited from the domain */
	typedef typename domain_t::scalar_t scalar_t;
	/** type of the local coordinate */
	typedef typename domain_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector */
	typedef Eigen::Matrix<scalar_t, num_nodes, 1> shape_t;
	/** \brief type of a \f$\nabla L(\xi)\f$ gradient matrix */
	typedef Eigen::Matrix<scalar_t, num_nodes, domain_t::dimension> dshape_t;
};

// Forward declaration
template <class Domain>
class constant_shape_set;

/** \brief Traits for constant shape sets */
template <class Domain>
struct shape_set_traits<constant_shape_set<Domain> >
{
	typedef Domain domain_t;	/**< \brief the domain type */
	static unsigned const num_nodes = 1;	/**< \brief number of nodes */
};

/**
* \brief Constant shape sets
*/
template <class Domain>
class constant_shape_set : public shape_set_base<constant_shape_set<Domain> >
{
public:
	typedef shape_set_base<constant_shape_set<Domain> > base_t;	/**< \brief the base class type */
	/** \brief the domain type */
	typedef typename base_t::domain_t domain_t;
	/** \brief type of the xi variable */
	typedef typename base_t::xi_t xi_t;
	/** \brief type of a shape function vector */
	typedef typename base_t::shape_t shape_t;
	/** \brief type of a shape function gradient matrix */
	typedef typename base_t::dshape_t dshape_t;

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

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return &(domain_t::get_center());
	}

	/** \brief return end iterator to the corner nodes
	* \return end iterator to corner nodes
	*/
	static xi_t const *corner_end(void)
	{
		return &(domain_t::get_center()) + 1;
	}

protected:
	static const shape_t m_shape;
	static const dshape_t m_dshape;
};


/** constant line shape set */
typedef constant_shape_set<line_domain> line_0_shape_set;
/** constant triangle shape set */
typedef constant_shape_set<tria_domain> tria_0_shape_set;
/** constant quad shape set */
typedef constant_shape_set<quad_domain> quad_0_shape_set;


// Forward declaration
template <class Domain>
class isoparam_shape_set;

/** \brief Traits for isoparametric shape sets */
template <class Domain>
struct shape_set_traits<isoparam_shape_set<Domain> >
{
	typedef Domain domain_t;	/**< \brief the domain type */
	static unsigned const num_nodes = domain_t::num_corners;	/**< \brief number of nodes */
};

/**
* \brief Isoparametric shape sets
*/
template <class Domain>
class isoparam_shape_set : public shape_set_base<isoparam_shape_set<Domain> >
{
public:
	typedef Domain domain_t;	/**< \brief the domain type */

	typedef shape_set_base<isoparam_shape_set<Domain> > base_t;	/**< \brief the base class type */
	/** \brief the xi location vector type */
	typedef typename base_t::xi_t xi_t;
	/** \brief type of a shape vector */
	typedef typename base_t::shape_t shape_t;
	/** \brief type of a shape gradient matrix */
	typedef typename base_t::dshape_t dshape_t;

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

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return domain_t::get_corners();
	}

	/** \brief return end iterator to the corner nodes
	* \return end iterator to corner nodes
	*/
	static xi_t const *corner_end(void)
	{
		return domain_t::get_corners() + domain_t::num_corners;
	}

	/**
	* \brief return corner at a given node number
	* \param [in] idx the node index
	* \return constant reference to the coordinates
	*/
	static xi_t const &corner_at(size_t idx)
	{
		return *(domain_t::get_corners() + idx);
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
	isoparam_shape_set<line_domain>::eval_dshape(typename isoparam_shape_set<line_domain>::xi_t const &)
{
	dshape_t dL;
	dL <<
		-0.5,
		+0.5;
	return dL;
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

// Forward declaration
class parallelogram_shape_set;

/** \brief Traits for parallelogram shapesets */
template<>
struct shape_set_traits<parallelogram_shape_set>
{
	typedef quad_domain domain_t;	/**< \brief the domain type */
	static unsigned const num_nodes = 3;	/**< \brief number of nodes */
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

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return domain_t::get_corners();
	}

	/** \brief return end iterator to the corner nodes
	* \return end iterator to corner nodes
	*/
	static xi_t const *corner_end(void)
	{
		return domain_t::get_corners() + domain_t::num_corners;
	}
};



// Forward declaration
class tria_2_shape_set;

/** \brief Traits for quadratic tria shapesets */
template<>
struct shape_set_traits<tria_2_shape_set>
{
	typedef tria_domain domain_t;	/**< \brief the domain type */
	static unsigned const num_nodes = 6;	/**< \brief number of nodes */
};

/**
* \brief quadratic 6-noded tria shape function set
*/
class tria_2_shape_set : public shape_set_base<tria_2_shape_set>
{
public:
	/**
	* \brief quadratic 6-noded tria shape functions
	* \param [in] xi the domain variable
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
	*/
	static dshape_t eval_dshape(xi_t const & _xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1];
		dshape_t dL;
		/** \todo these derivatives were computed by Matlab's simple(). Find an optimal evaluation for better performance. */
		dL <<
			4*eta+4*xi-3, 4*eta+4*xi-3,
			4-8*xi-4*eta, -4*xi,
			4*xi-1,       0,
			4*eta,        4*xi,
			0,            4*eta-1,
			-4*eta,       4-4*xi-8*eta;
		return dL;
	}

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return m_corners;
	}

	/** \brief return end iterator to the corner nodes
	* \return end iterator to corner nodes
	*/
	static xi_t const *corner_end(void)
	{
		return m_corners + num_nodes;
	}
    
    /** \todo this function should be defined generally, not for each new general shape set, because this is a sickness and we do not want it therefore */
    static xi_t const &corner_at(size_t idx)
	{
		return m_corners[idx];
	}
    

protected:
	static const xi_t m_corners[num_nodes];
};



// Forward declaration
class quad_2_shape_set;

/** \brief Traits for quadratic quad shapesets */
template<>
struct shape_set_traits<quad_2_shape_set>
{
	typedef quad_domain domain_t;	/**< \brief the domain type */
	static unsigned const num_nodes = 9;	/**< \brief number of nodes */
};

/**
* \brief quadratic 9-noded quad shape function set
*/
class quad_2_shape_set : public shape_set_base<quad_2_shape_set>
{
public:
	/**
	* \brief quadratic 9-noded quad shape functions
	* \param [in] xi the domain variable
	* \return the shape function vector
	*/
	static shape_t eval_shape(xi_t const &_xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1];
		shape_t L;
		L <<
			(1-xi)*( -xi)/2.0 * (1-eta)*( -eta)/2.0,
			(1-xi)*(1+xi)     * (1-eta)*( -eta)/2.0,
			(1+xi)*(  xi)/2.0 * (1-eta)*( -eta)/2.0,

			(1+xi)*(  xi)/2.0 * (1-eta)*(1+eta),

			(1+xi)*(  xi)/2.0 * (1+eta)*( +eta)/2.0,
			(1-xi)*(1+xi)     * (1+eta)*( +eta)/2.0,
			(1-xi)*( -xi)/2.0 * (1+eta)*( +eta)/2.0,

			(1-xi)*( -xi)/2.0 * (1-eta)*(1+eta),
			(1-xi)*(1+xi)     * (1-eta)*(1+eta);
		return L;
	}

	/**
	* \brief quadratic 9-noded quad shape function derivatives
	*/
	static dshape_t eval_dshape(xi_t const & _xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1], xi2 = xi*xi, eta2 = eta*eta;
		dshape_t dL;
		/** \todo these derivatives were computed by Matlab's simple(). Find an optimal evaluation for better performance. */
		dL <<
            eta*(2*xi-1)*(eta-1)/4, xi*(2*eta-1)*(xi-1)/4,
            -xi*eta*(eta-1),        -(xi2-1)*(2*eta-1)/2,
            eta*(2*xi+1)*(eta-1)/4, xi*(2*eta-1)*(xi+1)/4,
            -(2*xi+1)*(eta2-1)/2,   -xi*eta*(xi+1),
            eta*(2*xi+1)*(eta+1)/4, xi*(2*eta+1)*(xi+1)/4,
            -xi*eta*(eta+1),        -(xi2-1)*(2*eta+1)/2,
            eta*(2*xi-1)*(eta+1)/4, xi*(2*eta+1)*(xi-1)/4,
            -(2*xi-1)*(eta2-1)/2,   -xi*eta*(xi-1),
            2*xi*(eta2-1),          2*eta*(xi2-1);
		return dL;
	}

	/** \brief return begin iterator to the corner nodes
	* \return begin iterator to corner nodes
	*/
	static xi_t const *corner_begin(void)
	{
		return m_corners;
	}

	/** \brief return end iterator to the corner nodes
	* \return end iterator to corner nodes
	*/
	static xi_t const *corner_end(void)
	{
		return m_corners + num_nodes;
	}

	static xi_t const &corner_at(size_t idx)
	{
		return m_corners[idx];
	}

protected:
	static const xi_t m_corners[num_nodes];
};




template <>
typename constant_shape_set<quad_domain>::shape_t
	const constant_shape_set<quad_domain>::m_shape
	= constant_shape_set<quad_domain>::shape_t::Ones();

template <>
typename constant_shape_set<quad_domain>::dshape_t
	const constant_shape_set<quad_domain>::m_dshape
	= constant_shape_set<quad_domain>::dshape_t::Zero();

template <>
typename constant_shape_set<tria_domain>::shape_t
	const constant_shape_set<tria_domain>::m_shape
	= constant_shape_set<tria_domain>::shape_t::Ones();

template <>
typename constant_shape_set<tria_domain>::dshape_t
	const constant_shape_set<tria_domain>::m_dshape
	= constant_shape_set<tria_domain>::dshape_t::Zero();


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


tria_2_shape_set::xi_t
	const tria_2_shape_set::m_corners[tria_2_shape_set::num_nodes] = {
		tria_2_shape_set::xi_t(0.0, 0.0),
		tria_2_shape_set::xi_t(0.5, 0.0),
		tria_2_shape_set::xi_t(1.0, 0.0),
		tria_2_shape_set::xi_t(0.5, 0.5),
		tria_2_shape_set::xi_t(0.0, 1.0),
		tria_2_shape_set::xi_t(0.0, 0.5)
};


#endif // SHAPESET_HPP_INCLUDED

