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
		typedef int_ type;

		static int const value = N;
		static int const next = N+1;
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
	 * \brief metafunction returning difference of two integers
	 */
	template <int N, int M>
	struct less<int_<N>, int_<M> > : std::integral_constant<bool, N < M> {};
}

#endif

