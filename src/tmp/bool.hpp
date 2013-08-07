/**
 * \file bool.hpp
 * \ingroup tmp
 * \brief implementation of Boolean functions
 * \details The boolean type is alredy contained in type_traits as std::integral_constant<bool, x>.
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

	/** \brief disjunction of boolean constants */
	template <class...Args>
	struct or_ : std::false_type {};

	template <class...Args>
	struct or_<std::false_type, Args...> : or_<Args...> {};

	template <class...Args>
	struct or_<std::true_type, Args...> : std::true_type {};

	/** \brief conjunction of boolean constants */
	template <class...Args>
	struct and_ : std::false_type {};

	template <class...Args>
	struct and_<std::true_type, Args...> : and_<Args...> {};

	template <class...Args>
	struct and_<std::false_type, Args...> : std::false_type {};

	/**
	 * \brief IF control structure
	 * \tparam Choice a choice evaluated to a logical type
	 * \tparam T the type returned when Choice is true_type
	 * \tparam F the type returned when Choice is false_type
	 */
	template <class Cond, class T, class F>
	struct if_ : std::conditional<Cond::value, T, F> {};
}

#endif

