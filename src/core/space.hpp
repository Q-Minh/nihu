/** \file space.hpp
 * \ingroup funcspace
 * \brief declaration of class ::space
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu 
 */

#ifndef SPACE_HPP_INCLUDED
#define SPACE_HPP_INCLUDED

#include <Eigen/Dense>

/** \brief class representing a space with a scalar vaiable and a dimension
 * \tparam Scalar the scalar variable
 * \tparam Dimension the space dimensionality
 */
template <class Scalar, unsigned Dimension>
class space
{
public:
	/** \brief template parameter as nested type */
	typedef Scalar scalar_t;
	/** \brief template parameter as nested type */
	static unsigned const dimension = Dimension;

	/** \brief the location type */
	typedef Eigen::Matrix<scalar_t, Dimension, 1> location_t;
};

/** \brief specialisation for a 1D space of double */
typedef space<double, 1> space_1d;
/** \brief specialisation for a 2D space of double */
typedef space<double, 2> space_2d;
/** \brief specialisation for a 3D space of double */
typedef space<double, 3> space_3d;

#endif // SPACE_HPP_INCLUDED

