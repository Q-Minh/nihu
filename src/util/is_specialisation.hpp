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

#ifndef TYPE_INFO_HPP_INCLUDED
#define TYPE_INFO_HPP_INCLUDED

#include <type_traits>

namespace NiHu
{

template<class T>
struct enable_if_type { typedef void type; };

template<class T, class Enable = void>
struct is_specialisation : std::true_type {};

template<class T>
struct is_specialisation<T, typename enable_if_type<typename T::unspecialised>::type> : std::false_type
{};


}

#endif // TYPE_INFO_HPP_INCLUDED

