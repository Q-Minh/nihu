/** \file crtp_base.hpp
* \brief define CRTP helper functions
* \author Peter Fiala and Peter Rucz
* \ingroup util
*/

#ifndef CRTP_BASE_HPP_INCLUDED
#define CRTP_BASE_HPP_INCLUDED

/** \brief define CRTP helper function */
#define NIHU_CRTP_HELPERS \
	Derived const &derived() const { return static_cast<Derived const &>(*this); } \
	Derived &derived() { return static_cast<Derived &>(*this); }

#endif

