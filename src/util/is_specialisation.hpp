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
 * \file is_specialisation.hpp 
 * \brief Definition of the metafunction \ref is_specialisation
 * \ingroup util
 */

#ifndef IS_SPECIALISATION_HPP_INCLUDED
#define IS_SPECIALISATION_HPP_INCLUDED

#include <type_traits>

namespace NiHu
{

namespace internal
{
	/**
	 * \brief Helper metafunction for detecting specialisation 
	 * \details 
	 * The specialised case of the metafunction \ref is_specialisation tries to 
	 * instantiate this metafunction.
	 */
	template<class T>
	struct enable_if_type
	{ 
		typedef void type;
	};
} // end of namespace internal
	
/**
 * \brief Metafunction that determines if a type is a specialisation 
 * \returns True in the general case 
 * \details 
 * This is the general case, which is instantiated as a fallback when the 
 * specialised case below results in a substitution error,
 */
template<class T, class Enable = void>
struct is_specialisation
	: std::true_type {};


/**
 * \brief Metafunction that determines if a type is a specialisation 
 * \returns False in the case if T is not a specialisation
 * \details 
 * The specialised case tries to instantiate \ref internal::enable_if_type using
 * the type named as \c unspecialised in \c T. If the instantiation is 
 * successful, then \c T is not a specialisation.
 */
template<class T>
struct is_specialisation<T, 
	typename internal::enable_if_type<typename T::unspecialised>::type
>
	: std::false_type {};

} // end of namespace NiHu

#endif /* IS_SPECIALISATION_HPP_INCLUDED */
