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
	template <int N, int M>
	struct less<int_<N>, int_<M> > : std::integral_constant<bool, N < M> {};


	/** \brief compute maximum value of integral constants
	 * \todo this metafunction should be generalised for arbitrary types
	 */
	template <class Val, class...Args>
	struct max_
	{
		static unsigned const rest = max_<Args...>::value;
		static unsigned const x = Val::value;
		static unsigned const value = x > rest ? x : rest;
	};

	/** \brief specialisation of ::max_ for the one parameter case */
	template <class Val>
	struct max_<Val>
	{
		static unsigned const value = Val::value;
	};

}

#endif

