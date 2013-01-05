/**
 * \file algorithm.hpp
 * \brief vectoralgorithms
 */

#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "lambda.hpp"
#include "sequence.hpp"

/**
 * \brief accumulate elements of a container using a user-specified functor
 * \tparam Beg begin iterator
 * \tparam End end iterator
 * \tparam Inic initial value of accumulation
 * \tparam Fun accumulating functor, the default is plus
 */
template <class Beg, class End, class Init, class Fun = plus<_1,_2> >
struct accumulate : accumulate<
	typename next<Beg>::type,
	End,
	typename apply<Fun, Init, typename deref<Beg>::type>::type,
	Fun
> {};

/**
 * \brief terminating case of accumulate where the end iterator has been reached
 * \tparam End end iterator
 * \tparam Inic initial value of accumulation
 * \tparam Fun accumulating functor, the default is plus
 */
template <class End, class Inic, class Fun>
struct accumulate<End, End, Init, Fun> : Inic {};

#endif

