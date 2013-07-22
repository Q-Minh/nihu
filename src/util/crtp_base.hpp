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

