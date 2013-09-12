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
 * \file eigen_utils.hpp
 * \brief Collection of include files needed by Boonen13
 */

#ifndef EIGEN_UTILS_HPP_INCLUDED
#define EIGEN_UTILS_HPP_INCLUDED

#include <type_traits>
#include <Eigen/Dense>
#include <Eigen/StdVector>

/** \brief metafunction converting a type into an std::vector<T> type with Eigen allocator */
template <class T>
struct EigenStdVector
{
	typedef std::vector<T, Eigen::aligned_allocator<T> > type;
};

/** \brief metafunction determining if its argument is an Eigen expression or not
 * \tparam T the class to investigate
 */
template <class T>
struct is_eigen : std::is_base_of<
	Eigen::EigenBase<typename std::decay<T>::type>,
	typename std::decay<T>::type
> {};

template <class m1, class m2, int t>
struct is_eigen<Eigen::GeneralProduct<m1, m2, t> > : std::true_type {};

#endif // EIGEN_UTILS_HPP_INCLUDED

