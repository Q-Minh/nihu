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

/** \file potential_kernel.hpp
 * \ingroup library
 * \brief definition of tags describing various types of potentials
 */

#ifndef POTENTIAL_KERNEL_HPP_INCLUDED
#define POTENTIAL_KERNEL_HPP_INCLUDED

namespace potential
{
	/** \brief tag type to identify the Single Layer Potential */
	struct SLP { typedef SLP type; };
	/** \brief tag type to identify the Double Layer Potential */
	struct DLP { typedef DLP type; };
	/** \brief tag type to identify the transposed Double Layer Potential */
	struct DLPt { typedef DLPt type; };
	/** \brief tag type to identify the Hypersingular Potential */
	struct HSP { typedef HSP type; };
}

#endif // POTENTIAL_KERNEL_HPP_INCLUDED

