/**
 * \file lambda.hpp
 * \brief implementation of placeholders and lambda functions
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef LAMBDA_HPP_INCLUDED
#define LAMBDA_HPP_INCLUDED

#include "bool.hpp"

namespace tmp
{
	/**
	 * \brief placeholder that selects N-th argument
	 * \tparam N argument index
	 */
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


	/**
	 * \brief metafunction returning std::true_type if its argument is a placeholder
	 */
	template <class C>
	struct isPlaceholder : std::false_type {};

	template <int N>
	struct isPlaceholder<arg<N> > : std::true_type {};


	/**
	 * \brief metafunction returning std::true_type if its argument is a placeholder expression
	 */
	template <class C>
	struct isPlaceholderExpression : std::false_type {};

	template <int N>
	struct isPlaceholderExpression<arg<N> > : std::true_type {};

	template <template <class Arg1> class MF, class Arg1>
	struct isPlaceholderExpression<MF<Arg1> > : isPlaceholderExpression<Arg1> {};

	template <template <class Arg1, class Arg2> class MF, class Arg1, class Arg2>
	struct isPlaceholderExpression<MF<Arg1, Arg2> > : or_<
		typename isPlaceholderExpression<Arg1>::type,
		typename isPlaceholderExpression<Arg2>::type
	> {};

	template <template <class Arg1, class Arg2, class Arg3> class MF, class Arg1, class Arg2, class Arg3>
	struct isPlaceholderExpression<MF<Arg1, Arg2, Arg3> > : or_<
		typename isPlaceholderExpression<Arg1>::type,
		typename isPlaceholderExpression<Arg2>::type,
		typename isPlaceholderExpression<Arg3>::type
	> {};


	namespace internal
	{
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
	}

	/**
	 * \brief generate metafunction class from placeholder expression
	 * \tparam Exp placeholder expression or metafunction class
	 * \return metafunction class encapsulating the placeholder expression or the input class itself
	 */
	template <class Exp>
	struct lambda : if_<
		typename isPlaceholderExpression<Exp>::type,
		typename internal::lambda_plExp<Exp>::type,
		Exp
	> {};

	template <class Fun, class Arg1 = void_, class Arg2 = void_>
	struct apply : lambda<Fun>::type::template apply<Arg1, Arg2> {};
}

#endif

