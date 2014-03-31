// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
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

/** \file blind_singular_quadrature.hpp
 * \ingroup quadrature
 */

#ifndef BLIND_SINGULAR_QUADRATURE_HPP_INCLUDED
#define BLIND_SINGULAR_QUADRATURE_HPP_INCLUDED

#include "formalism.hpp"
#include "blind_transform_selector.hpp"
#include "duffy_quadrature.hpp"

template <class BlindTransform, class RegularFamily, class LSet>
struct blind_singular_quadrature;

template <class RegularFamily, class LSet>
struct blind_singular_quadrature<blind_transform::duffy, RegularFamily, LSet>
{
	typedef duffy_quadrature<RegularFamily, LSet> type;
};

#endif // BLIND_SINGULAR_QUADRATURE_HPP_INCLUDED
