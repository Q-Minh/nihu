/**
 * \file bool.hpp
 * \ingroup tmp
 * \brief implementation of Boolean functions
 * \details The boolean type is alredy contained in <type_traits> as std::integral_constant<bool, x>.
 * This file implements some Boolean functions and the compile time if_ control structure
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef BOOL_HPP_INCLUDED
#define BOOL_HPP_INCLUDED

#include <type_traits>

namespace tmp
{
	/**
	 * \brief Boolean negation
	 * \tparam A Boolean data
	 * \returns not A
	 */
	template <class A>
	struct not_ : std::integral_constant<bool, !A::value> {};

	/** \brief Boolean disjunction of up to three Boolean data */
	template <class A, class B, class C = std::false_type>
	struct or_ : std::integral_constant<bool, A::value || B::value || C::value> {};

	/** \brief Boolean conjunction of up to three Boolean data */
	template <class A, class B, class C = std::true_type>
	struct and_ : std::integral_constant<bool, A::value && B::value && C::value> {};

	/**
	 * \brief IF control structure
	 * \tparam Choice a choice evaluated to a logical type
	 * \tparam T the type returned when Choice is true_type
	 * \tparam F the type returned when Choice is false_type
	 */
	template <class Cond, class T, class F>
	struct if_;

	/** \brief specialisation of if_ for the true choice case */
	template <class T, class F>
	struct if_<std::true_type, T, F> { typedef T type; };

	/** \brief specialisation of if_ for the false choice case */
	template <class T, class F>
	struct if_<std::false_type, T, F> { typedef F type; };
}

#endif

