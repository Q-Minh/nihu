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

/** \file crtp_base.hpp
* \brief define CRTP helper functions and metafunctions
* \author Peter Fiala and Peter Rucz
* \ingroup util
*/

#ifndef CRTP_BASE_HPP_INCLUDED
#define CRTP_BASE_HPP_INCLUDED

/** \brief define CRTP helper function */
#define NIHU_CRTP_HELPERS \
	Derived const &derived() const { return static_cast<Derived const &>(*this); } \
	Derived &derived() { return static_cast<Derived &>(*this); }


/** \brief metafunction returning its first argument and ignoring all subsequent
 * \details used for CRTP decltype
 */
template <class T, class...Ignore>
struct ignore
{
	typedef T type;
};

/** \brief crtp decltype helper function */
template <class Derived, class Dummy>
typename ignore<Derived, Dummy>::type const* const_crtp_ptr(void)
{
	return static_cast<typename ignore<Derived, Dummy>::type const*>(nullptr);
}

#endif // CRTP_BASE_HPP_INCLUDED

