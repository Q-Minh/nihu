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

/**
 * \file global_definitions.hpp
 * \brief global constants governing some accuracy parameters
 * \ingroup core
 */

#ifndef GLOBAL_DEFINITIONS_HPP_INCLUDED
#define GLOBAL_DEFINITIONS_HPP_INCLUDED

#include "field_type_acceleration_option.hpp"

/** \brief acceleration::soft or acceleration::hard */
typedef acceleration::hard GLOBAL_ACCELERATION;
/**
 * \brief the maximal order of accelerated quadratures and field type accelerators
 * \todo increase Dunavant order to avoid overindexing
 */
unsigned const GLOBAL_MAX_ORDER = 9;
/** \brief the global accuracy of integraions */
unsigned const GLOBAL_ACCURACY = 3;

#endif

