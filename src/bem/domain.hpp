/**
 * \file domain.hpp
 * \brief declaration of class domain and its derived domains
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include <Eigen/Dense>

/**
 * \brief class Domain represents a subset of the xi space. All elements are defined on a domain.
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
};

/** \brief a 1D line domain \f$-1 \le \xi \le +1\f$*/
class line_domain : public domain<1, 2> {};
/** \brief a 2D triangle domain */
class tria_domain : public domain<2, 3> {};
/** \brief a 2D quad domain */
class quad_domain : public domain<2, 4> {};
/** \brief a 3D brick domain */
class brick_domain : public domain<3, 8> {};

#endif

