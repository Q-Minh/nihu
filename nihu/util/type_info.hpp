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
 * \file type_info.hpp 
 * \brief Utility functions for attaining type informations
 * \ingroup util
 */

#ifndef TYPE_INFO_HPP_INCLUDED
#define TYPE_INFO_HPP_INCLUDED

#include <iostream>
#include <type_traits>

namespace NiHu
{

/**
 * \brief Print type information 
 * \tparam T The requested type 
 * \param os Output stream to print to, default is \c std::cout
 * \return The output stream
 */
template <class T>
std::ostream & print_type_info(std::ostream &os = std::cout)
{
	os << (std::is_const<typename std::remove_reference<T>::type >::value ? "const " : "");
	os << (std::is_same<typename std::decay<T>::type, int>::value ? "int " : "");
	os << (std::is_same<typename std::decay<T>::type, char>::value ? "char " : "");
	os << (std::is_lvalue_reference<T>::value ? "&" : "");
	os << (std::is_rvalue_reference<T>::value ? "&&" : "");
	return os;
}

} // end of namespace NiHu

#endif /* TYPE_INFO_HPP_INCLUDED */
