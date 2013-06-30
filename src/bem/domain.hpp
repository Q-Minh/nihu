/**
 * \file domain.hpp
 * \brief declaration of class ::domain and its derived domains
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include "space.hpp"

/** \brief metafunction assigning an id to a domain */
template <class domain_t>
struct domain_id;


/**
 * \brief a subset of the \f$\xi\f$ space. All elements are defined on a domain.
 * \tparam Scalar the scalar type of the space
 * \tparam Dimension the dimensionality of the space
 * \tparam NumCorners the number of corners of the domain
 */
template <class Space, unsigned NumCorners>
class domain
{
public:
	/** \brief template parameter as nested type */
	typedef Space space_t;
	/** \brief template argument as nested type */
	static unsigned const num_corners = NumCorners;
	
	/** \brief scalar type inherited from the space */
	typedef typename space_t::scalar_t scalar_t;
	/** \brief dimension inherited from the sapce */
	static unsigned const dimension = space_t::dimension;

	/** \brief the domain id */
	static unsigned const id = domain_id<domain>::value;

	/** \brief location vector renamed */
	typedef typename space_t::location_t xi_t;

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


/** \brief metafunction assigning an id to a domain */
template <class domain_t>
struct domain_id
{
	/** \brief default computed value of a domain id */
	static unsigned const value =
		domain_t::dimension * 10 +
		domain_t::num_corners;
};



/** \brief a 1D line domain \f$-1 \le \xi \le +1\f$*/
typedef domain<space_1d, 2> line_domain;

/** \brief a 2D triangle domain */
typedef domain<space_2d, 3> tria_domain;

/** \brief a 2D quad domain */
typedef domain<space_2d, 4> quad_domain;

/** \brief a 3D brick domain */
typedef domain<space_3d, 8> brick_domain;


template<>
line_domain::xi_t
	const line_domain::m_center =
	line_domain::xi_t::Zero();


template<>
line_domain::corners_t
	const line_domain::m_corners = {
	line_domain::xi_t::Constant(-1.0),
	line_domain::xi_t::Constant(1.0)
	};


template<>
tria_domain::xi_t
	const tria_domain::m_center =
	tria_domain::xi_t::Constant(1.0/3.0);


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

