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

/// \file mex_matrix.hpp
/// \brief A Matlab mex matrix interface
/// \ingroup matlab
/// \details 
/// The interface makes it possible to use Matlab-borne matrices in C++ and to 
/// create Matlab matrices in C++. The interface handles real and complex
/// matrices in a convenient manner, hiding mex implementation details from the
/// C++ programmer.
 
#ifndef MEX_MATRIX_HPP_INCLUDED
#define MEX_MATRIX_HPP_INCLUDED

#include <mex.h>
#include <matrix.h>

#if MX_HAS_INTERLEAVED_COMPLEX
#include "mex_matrix_interleaved.hpp"
#else
#include "mex_matrix_separate.hpp"
#endif

#endif /// MEX_MATRIX_HPP_INCLUDED
