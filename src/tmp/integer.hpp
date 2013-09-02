/**
 * \file integer.hpp
 * \brief integer type representation and basic integer arithmetics
 */

#ifndef INTEGER_HPP
#define INTEGER_HPP

#include "operator.hpp"
#include "bool.hpp"

namespace tmp
{
	/**
	 * \brief integer type representation
	 */
	template <int N>
	struct int_
	{
		/** \brief self returning metafunction */
		typedef int_ type;
		/** \brief the stored value */
		static int const value = N;
		/** \brief the next value */
		static int const next = N+1;
		/** \brief the previous value */
		static int const prev = N-1;
	};

	/**
	 * \brief metafunction returning next integer
	 */
	template <int N>
	struct next<int_<N> > : int_<int_<N>::next> {};

	/**
	 * \brief metafunction returning previous integer
	 */
	template <int N>
	struct prev<int_<N> > : int_<int_<N>::prev> {};

	/**
	 * \brief metafunction returning sum of two integers
	 */
	template <int N, int M>
	struct plus<int_<N>, int_<M> > : int_<N+M> {};

	/**
	 * \brief metafunction returning difference of two integers
	 */
	template <int N, int M>
	struct minus<int_<N>, int_<M> > : int_<N-M> {};

	/**
	 * \brief metafunction returning difference of two integers
	 */
	template <int N, int M>
	struct mul<int_<N>, int_<M> > : int_<N*M> {};

	/**
	 * \brief metafunction returning true if first is less than second
	 */
	template <class N, class M> struct less;

	template <int N, int M>
	struct less<int_<N>, int_<M> > : std::integral_constant<bool, N < M> {};

	/**
	 * \brief metafunction returning true if first is greater than second
	 */
	template <class N, class M> struct greater;

	template <int N, int M>
	struct greater<int_<N>, int_<M> > : std::integral_constant<bool, !(N <= M)> {};


	/** \brief compute maximum value of integral constants
	 * \todo should be generalised for general type
	 */
	template <class Val, class...Args>
	struct max_ : max_<Val, typename max_<Args...>::type> {};

	/** \brief specialisation of max_ for the two parameter case */
	template <class Val1, class Val2>
	struct max_<Val1, Val2>
	{
		/** \brief the maximum value */
		static int const value = Val1::value > Val2::value ? Val1::value : Val2::value;
		/** \brief the maximum type */
		typedef int_<value> type;
	};

	/** \brief compute minimum value of integral constants
	 * \todo should be generalised for general type
	 */
	template <class Val, class...Args>
	struct min_ : min_<Val, typename min_<Args...>::type> {};

	/** \brief specialisation of min_ for the two parameter case */
	template <class Val1, class Val2>
	struct min_<Val1, Val2>
	{
		/** \brief the minimum value */
		static int const value = Val1::value < Val2::value ? Val1::value : Val2::value;
		/** \brief the minimum type */
		typedef int_<value> type;
	};
}

#endif

