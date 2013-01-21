/**
 * \file bool.hpp
 * \brief implementation of bool_ and Boolean functions
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef BOOL_HPP_INCLUDED
#define BOOL_HPP_INCLUDED

namespace tmp
{
	/**
	 * \brief Boolean type representation
	 * \tparam x the underlying Boolean value
	 */
	template <bool x>
	struct bool_
	{
		static bool const value = x;
		typedef bool_ type;
	};

	/** \brief true type */
	typedef bool_<true> true_;
	/** \brief false type */
	typedef bool_<false> false_;

	template <class Cond, class T, class F>
	struct if_;

	template <class T, class F>
	struct if_<true_, T, F> { typedef T type; };

	template <class T, class F>
	struct if_<false_, T, F> { typedef F type; };

	template <class A>
	struct not_ : bool_<!A::value> {};
	
	template <class A, class B, class C = false_>
	struct or_ : bool_<A::value || B::value || C::value> {};
	
	template <class A, class B, class C = true_>
	struct and_ : bool_<A::value && B::value && C::value> {};
}

#endif

