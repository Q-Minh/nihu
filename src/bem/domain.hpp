/**
 * \file domain.hpp
 * \brief declaration of class domain and its derived domains
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \details Domains are the domain set of geometrical transformations that desrcibe elements.
 * Each element is mapped from its domain set to its actual geometry by means
 * of a transformation \f$x = x(\xi)\f$, where \f$\xi \in \mathcal{D}\f$. In this description,
 * \f$\mathcal{D}\f$ denotes the element's domain.
 */
#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include <Eigen/Dense>

/**
 * \brief a subset of the \f$\xi\f$ space. All elements are defined on a domain.
 * \tparam Dimension the dimensionality of the space
 * \tparam ID A domain ID used when elements are inserted into a mesh
 */
template <unsigned Dimension, unsigned ID>
class domain
{
public:
	/** \brief template argument as nested type */
	static unsigned const dimension = Dimension;
	/** \brief template argument as nested type */
	static unsigned const id = ID;

	/** \brief the scalar type of the location vector */
    typedef double scalar_t;
	/** \brief the location vector */
	typedef Eigen::Matrix<double, dimension, 1> xi_t;

	/**
	 * \brief return the central point of the domain
	 * \return center point
	 */
	static xi_t const &get_center(void)
	{
		return m_center;
	}

protected:
	/** \brief the center point of the domain */
	static xi_t m_center;
};

/** \brief a 1D line domain \f$-1 \le \xi \le +1\f$*/
typedef domain<1, 2> line_domain;
template<>
line_domain::xi_t line_domain::m_center = line_domain::xi_t::Zero();

/** \brief a 2D triangle domain */
typedef domain<2, 3> tria_domain;
template<>
tria_domain::xi_t tria_domain::m_center = tria_domain::xi_t::Ones()/3.0;

/** \brief a 2D quad domain */
typedef domain<2, 4> quad_domain;
template<>
quad_domain::xi_t quad_domain::m_center = quad_domain::xi_t::Zero();

/** \brief a 3D brick domain */
typedef domain<3, 8> brick_domain;
template<>
brick_domain::xi_t brick_domain::m_center = brick_domain::xi_t::Zero();

#endif

