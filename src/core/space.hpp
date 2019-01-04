// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2019  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2019  Peter Rucz <rucz@hit.bme.hu>
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

namespace NiHu
{

/** \brief class representing a coordinate space with a scalar and a dimension
 * \tparam Scalar the scalar variable
 * \tparam Dimension the space dimensionality
 */
template <class Scalar, unsigned Dimension>
struct space
{
    /** self-returning */
    typedef space type;

	/** \brief template parameter as nested type */
	typedef Scalar scalar_t;

	/** \brief integral constants */
	enum {
		/** \brief template parameter as nested type */
		dimension = Dimension
	};

	/** \brief the location type */
	typedef Eigen::Matrix<scalar_t, Dimension, 1> location_t;
};

/** \brief specialisation for a 1D space */
template <class Scalar = double>
using space_1d = space<Scalar, 1>;

/** \brief specialisation for a 2D space */
template <class Scalar = double>
using space_2d = space<Scalar, 2>;

/** \brief specialisation for a 3D space */
template <class Scalar = double>
using space_3d = space<Scalar, 3>;

} // end of namespace NiHu

#endif // SPACE_HPP_INCLUDED

