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

/**
 * \file field_type_acceleration_option.hpp
 * \ingroup quadrature
 * \brief definition of field type acceleration options
 */

#ifndef FIELD_TYPE_ACCELERATION_OPTIONS_HPP_INCLUDED
#define FIELD_TYPE_ACCELERATION_OPTIONS_HPP_INCLUDED

/** \brief collection of acceleration-types */
namespace acceleration
{
	/** \brief real acceleration */
	struct hard {};
	/** \brief view-acceleration */
	struct soft{};
}

#endif // FIELD_TYPE_ACCELERATION_OPTIONS_HPP_INCLUDED
