/**
 * \file bool.hpp
 * \brief implementation of bool_ and Boolean functions
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef BOOL_HPP_INCLUDED
#define BOOL_HPP_INCLUDED

#include <type_traits>

namespace tmp
{
	template <class Cond, class T, class F>
	struct if_;

	template <class T, class F>
	struct if_<std::true_type, T, F> { typedef T type; };

	template <class T, class F>
	struct if_<std::false_type, T, F> { typedef F type; };

	template <class A>
	struct not_ : std::integral_constant<bool, !A::value> {};

	template <class A, class B, class C = std::false_type>
	struct or_ : std::integral_constant<bool, A::value || B::value || C::value> {};

	template <class A, class B, class C = std::true_type>
	struct and_ : std::integral_constant<bool, A::value && B::value && C::value> {};
}

#endif

