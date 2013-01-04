#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "lambda.hpp"
#include "sequence.hpp"

template <class Beg, class End, class Inic, class Fun = plus<_1,_2> >
struct accumulate : accumulate<
	typename next<Beg>::type,
	End,
	typename apply<Fun, Inic, typename deref<Beg>::type>::type,
	Fun
> {};

template <class End, class Inic, class Fun>
struct accumulate<End, End, Inic, Fun> : Inic {};

#endif

