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

/** \file space.hpp
 * \ingroup funcspace
 * \brief declaration of class ::space
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

