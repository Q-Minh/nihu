/**
 * \file domain.hpp
 * \brief declaration of class ::domain and its derived domains
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \details Domains are the domain set of geometrical transformations that desrcibe elements.
 * Each element is mapped from its domain set to its actual geometry by means
 * of a transformation \f$x = x(\xi)\f$, where \f$\xi \in \mathcal{D}\f$. In this description,
 * \f$\mathcal{D}\f$ denotes the element's domain.
 */
#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include "includes.h"

/**
 * \brief a subset of the \f$\xi\f$ space. All elements are defined on a domain.
 * \tparam Dimension the dimensionality of the space
 * \tparam NumCorners the number of corners of the domain
 */
template <unsigned Dimension, unsigned NumCorners>
class domain
{
public:
	/** \brief template argument as nested type */
	static unsigned const dimension = Dimension;
	/** \brief template argument as nested type */
	static unsigned const num_corners = NumCorners;
	/** \todo in the present implementation id equals the number of corner nodes. fix this! */
	static unsigned const id = num_corners;

	/** \brief scalar type of the location vector */
    typedef double scalar_t;
	/** \brief location vector */
	typedef Eigen::Matrix<scalar_t, dimension, 1> xi_t;
	/** \brief type of the corners' array */
	typedef const xi_t corners_t[num_corners];

	/**
	 * \brief return the central point of the domain
	 * \return center point
	 */
	static xi_t const &get_center(void)
	{
		return m_center;
	}

	/**
	 * \brief return pointer to array of corners
	 * \return pointer to array of corners
	 */
	static xi_t const *get_corners(void)
	{
		return m_corners;
	}
	
	/**
	 * \brief return reference to a coner point
	 * \return constant reference to a corner point
	 */
	static xi_t const &get_corner(unsigned idx)
	{
		return get_corners()[idx];
	}

protected:
	/** \brief the center point of the domain */
	static xi_t const m_center;
	/** \brief the corner points of the domain */
	static corners_t const m_corners;
};

/** \brief a 1D line domain \f$-1 \le \xi \le +1\f$*/
typedef domain<1, 2> line_domain;

/** \brief a 2D triangle domain */
typedef domain<2, 3> tria_domain;

/** \brief a 2D quad domain */
typedef domain<2, 4> quad_domain;

/** \brief a 3D brick domain */
typedef domain<3, 8> brick_domain;






template<>
line_domain::xi_t
	const line_domain::m_center =
	line_domain::xi_t::Zero();


template<>
line_domain::corners_t
	const line_domain::m_corners = {
	-line_domain::xi_t::Ones(),
	line_domain::xi_t::Ones()
	};


template<>
tria_domain::xi_t
	const tria_domain::m_center =
	tria_domain::xi_t::Ones()/3.0;


template<>
tria_domain::corners_t const tria_domain::m_corners = {
	tria_domain::xi_t(0.0,0.0),
	tria_domain::xi_t(1.0,0.0),
	tria_domain::xi_t(0.0,1.0)
	};


template<>
quad_domain::xi_t
	const quad_domain::m_center =
	quad_domain::xi_t::Zero();


template<>
quad_domain::corners_t const quad_domain::m_corners = {
	quad_domain::xi_t(-1.0,-1.0),
	quad_domain::xi_t( 1.0,-1.0),
	quad_domain::xi_t( 1.0, 1.0),
	quad_domain::xi_t(-1.0, 1.0)
	};


template<>
brick_domain::corners_t
	const brick_domain::m_corners = {
	brick_domain::xi_t(-1.0,-1.0,-1.0),
	brick_domain::xi_t( 1.0,-1.0,-1.0),
	brick_domain::xi_t( 1.0, 1.0,-1.0),
	brick_domain::xi_t(-1.0, 1.0,-1.0),
	brick_domain::xi_t(-1.0,-1.0, 1.0),
	brick_domain::xi_t( 1.0,-1.0, 1.0),
	brick_domain::xi_t( 1.0, 1.0, 1.0),
	brick_domain::xi_t(-1.0, 1.0, 1.0)
};


template<>
brick_domain::xi_t
	const brick_domain::m_center =
	brick_domain::xi_t::Zero();




#endif

