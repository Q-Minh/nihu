/**
 * \file lambda.hpp
 * \brief lambda function implementation
 */

#ifndef LAMBDA_HPP
#define LAMBDA_HPP

#include "placeholder.hpp"

template <class Fun>
struct lambda_plExp // dummy case
{
	typedef struct {
		template <class A1 = void_, class A2 = void_>
		struct apply { typedef void_ type; }; // GCC needs this to compile, but this class is never used
	} type;
};

template <int N>
struct lambda_plExp<arg<N> >
{
	typedef struct
	{
		template <class A1 = void_, class A2 = void_>
		struct apply : arg<N>::template apply<A1, A2> {};
	} type;
};

template <template <class a> class MetaFun, class a1>
struct lambda_plExp<MetaFun<a1> >
{
	typedef struct
	{
		template <class A1 = void_, class A2 = void_>
		struct apply : MetaFun< typename lambda_plExp<a1>::type::template apply<A1, A2>::type > {};
	} type;
};

template <template <class a, class b> class MetaFun, class a1, class a2>
struct lambda_plExp<MetaFun<a1, a2> >
{
	typedef struct
	{
		template <class A1 = void_, class A2 = void_>
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

template <template <class a1, class a2, class a3> class MetaFun, class a1, class a2, class a3>
struct lambda_plExp<MetaFun<a1, a2, a3> >
{
	typedef struct
	{
		template <class A1 = void_, class A2 = void_>
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
			>::type,
			typename if_<
				typename isPlaceholderExpression<a3>::type,
				typename lambda_plExp<a3>::type::template apply<A1, A2>::type,
				a3
			>::type
		> {};
	} type;
};

/**
 * \brief generate metafunction class from placeholder expression
 * \tparam Exp placeholder expression or metafunction class
 * \return metafunction class encapsulating the placeholder expression or the input class itself
 */
template <class Exp>
struct lambda : if_<
	typename isPlaceholderExpression<Exp>::type,
	typename lambda_plExp<Exp>::type,
	Exp
> {};

template <class Fun, class Arg1 = void_, class Arg2 = void_>
struct apply : lambda<Fun>::type::template apply<Arg1, Arg2> {};

#endif

