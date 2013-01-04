#ifndef PLACEHOLDER_HPP
#define PLACEHOLDER_HPP

#include "bool.hpp"

template <int N>
struct arg;

struct void_;

template <>
struct arg<1>
{
	template <class A1 = void_, class A2 = void_>
	struct apply
	{
		typedef A1 type;
	};
};

typedef arg<1> _1;

template <>
struct arg<2>
{
	template <class A1 = void_, class A2 = void_>
	struct apply
	{
		typedef A2 type;
	};
};

typedef arg<2> _2;


template <class C>
struct isPlaceholder : false_ {};

template <int N>
struct isPlaceholder<arg<N> > : true_ {};


template <class C>
struct isPlaceholderExpression : false_ {};

template <int N>
struct isPlaceholderExpression<arg<N> > : true_ {};

template <template <class Arg1> class MF, class Arg1>
struct isPlaceholderExpression<MF<Arg1> > : isPlaceholderExpression<Arg1> {};

template <template <class Arg1, class Arg2> class MF, class Arg1, class Arg2>
struct isPlaceholderExpression<MF<Arg1, Arg2> > : or_<typename isPlaceholderExpression<Arg1>::type, typename isPlaceholderExpression<Arg2>::type> {};

#endif
