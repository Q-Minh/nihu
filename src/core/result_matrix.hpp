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
 * \file result_matrix.hpp
 * \brief Encapsulates different result matrix types
 * \ingroup assembly
 */

#ifndef RESULT_MATRIX_HPP_INCLUDED
#define RESULT_MATRIX_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "../util/eigen_utils.hpp"

namespace NiHu
{

template <class T>
struct is_result_matrix_impl :
	std::false_type {};

template <class T>
struct is_result_matrix_impl<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > :
	std::true_type {};

template <class Mat>
struct is_result_matrix : is_result_matrix_impl<
	typename std::decay<Mat>::type
> {};

} // end of namespace NiHu

#endif // RESULT_MATRIX_HPP_INCLUDED
