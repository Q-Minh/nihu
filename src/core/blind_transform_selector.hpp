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

/** \file blind_transform_selector.hpp
 * \brief select a blind transform method to a singularity type and a reference domain
 * \ingroup quadrature
 */

#ifndef BLIND_TRANSFORM_SELECTOR_HPP_INCLUDED
#define BLIND_TRANSFORM_SELECTOR_HPP_INCLUDED

#include "asymptotic_types.hpp"
#include "domain.hpp"
#include "library/lib_domain.hpp"

/** \brief defining blind transform algorithms */
namespace blind_transform
{
	/** \brief Duffy type polar transformation */
	struct duffy {};
	/** \brief t^2 transformation to cancel log singularities */
	struct square {};
}

/** \brief assign a blind transformation method to a singularity type and a renference domain
 * \tparam asymptotic_type the singularity type
 * \tparam domain the reference domain
 */
template <class asymptotic_type, class domain>
struct blind_transform_selector;

/** \brief assign a blind transformation method to log<r> singularity
 */
template <>
struct blind_transform_selector<asymptotic::log<1>, line_domain>
{
	typedef blind_transform::square type;
};

/** \brief assign a blind transformation method to 1/r singularity and tria domain
 */
template <>
struct blind_transform_selector<asymptotic::inverse<1>, tria_domain>
{
	typedef blind_transform::duffy type;
};

/** \brief assign a blind transformation method to 1/r singularity and quad domain
 */
template <>
struct blind_transform_selector<asymptotic::inverse<1>, quad_domain>
{
	typedef blind_transform::duffy type;
};


#endif // BLIND_TRANSFORM_SELECTOR_INCLUDED

