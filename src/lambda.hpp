#ifndef LAMBDA_HPP
#define LAMBDA_HPP

#include "placeholder.hpp"

template <class Fun>
struct lambda_plExp // dummy case
{
	typedef struct {
		template <class A1, class A2>
		struct apply; // intentionally left empty, never instantiated
	} type;
};

template <int N>
struct lambda_plExp<arg<N> >
{
	typedef struct
	{
		template <class A1, class A2>
		struct apply : arg<N>::template apply<A1, A2> {};
	} type;
};

template <template <class a> class MetaFun, class a1>
struct lambda_plExp<MetaFun<a1> >
{
	typedef struct
	{
		template <class A1, class A2>
		struct apply : MetaFun< typename lambda_plExp<a1>::type::template apply<A1, A2>::type > {};
	} type;
};

template <template <class a, class b> class MetaFun, class a1, class a2>
struct lambda_plExp<MetaFun<a1, a2> >
{
	typedef struct
	{
		template <class A1, class A2>
		struct apply : MetaFun<
			typename if_<
				typename isPlaceholderExpression<a1>::type,
				typename lambda_plExp<a1>::type::template apply<A1, A2>::type,
				a1
			>::type,
			typename if_<
				typename isPlaceholderExpression<a2>::type,
				typename lambda_plExp<a2>::type::template apply<A1, A2>::type,
				a2
			>::type
		> {};
	} type;
};

template <class Fun>
struct lambda : if_<
	typename isPlaceholderExpression<Fun>::type,
	typename lambda_plExp<Fun>::type,
	Fun
> {};

#endif
