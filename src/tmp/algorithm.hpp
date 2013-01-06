/**
 * \file algorithm.hpp
 * \brief vectoralgorithms
 */

#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "lambda.hpp"
#include "sequence.hpp"
#include "operator.hpp"

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
 * \tparam Init initial value of accumulation
 * \tparam Fun accumulating functor
 */
template <class End, class Init, class Fun>
struct accumulate<End, End, Init, Fun> : Init {};

template <class Beg, class End>
struct min : accumulate<Beg, End, typename deref<Beg>::type, if_< less<_1,_2>, _1, _2> > {};

template <class Beg, class End>
struct max : accumulate<Beg, End, typename deref<Beg>::type, if_< less<_1,_2>, _2, _1> > {};

#endif

