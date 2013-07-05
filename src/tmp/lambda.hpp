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
	template <unsigned N, class T, class ...Args>
	struct select_argument : select_argument<N-1, Args...> {};

	template <class T, class ...Args>
	struct select_argument<1U, T, Args...>
	{
		typedef T type;
	};

	/**
	 * \brief placeholder that selects N-th argument
	 * \tparam N argument index
	 */
	template <unsigned N>
	struct arg
	{
		template <class...Args>
		struct apply : select_argument<N, Args...> {};
	};
	
	typedef arg<1> _1;
	typedef arg<2> _2;
	typedef arg<3> _3;
	

	/** \brief used to represent no parameter */
	struct void_;

	/**
	 * \brief metafunction returning std::true_type if its argument is a placeholder
	 */
	template <class C>
	struct isPlaceholder : std::false_type {};

	template <unsigned N>
	struct isPlaceholder<arg<N> > : std::true_type {};


	/**
	 * \brief metafunction returning std::true_type if its argument is a placeholder expression
	 */
	template <class C>
	struct isPlaceholderExpression : std::false_type {};

	template <unsigned N>
	struct isPlaceholderExpression<arg<N> > : std::true_type {};

	template <class T, class...Args>
	struct containsPlaceholderExpression;
	
	template <class First, class...Args>
	struct containsPlaceholderExpressionImpl;
	 
	template <class...Args>
	struct containsPlaceholderExpressionImpl<std::true_type, Args...> : std::true_type {};
	 
	template <>
	struct containsPlaceholderExpressionImpl<std::false_type> : std::false_type {};
	 
	template <class...Args>
	struct containsPlaceholderExpressionImpl<std::false_type, Args...> : containsPlaceholderExpression<Args...> {};
	 
	template <class T, class...Args>
	struct containsPlaceholderExpression : containsPlaceholderExpressionImpl<
		typename isPlaceholderExpression<T>::type,
		Args...
	> {};

	template <template <class...Args> class MF, class...Args>
	struct isPlaceholderExpression<MF<Args...> > :
		containsPlaceholderExpression<Args...> {};


	namespace internal
	{
		template <class Fun>
		struct lambda_plExp // dummy case
		{
			typedef struct {
				template <class A1 = void_, class A2 = void_>
				// GCC needs this to compile, but this class is never used
				struct apply { typedef void_ type; };
			} type;
		};

		template <unsigned N>
		struct lambda_plExp<arg<N> >
		{
			typedef struct
			{
				template <class...Args>
				struct apply : arg<N>::template apply<Args...> {};
			} type;
		};

		template <template <class a> class MetaFun, class a1>
		struct lambda_plExp<MetaFun<a1> >
		{
			typedef struct
			{
				template <class...Args>
				struct apply
					: MetaFun<
						typename lambda_plExp<a1>::type::template apply<Args...>::type
						>
				{};
			} type;
		};

		template <class BOOL, class a, class...Args>
		struct cond_eval_arg : lambda_plExp<a>::type::template apply<Args...> {};

		template <class a, class...Args>
		struct cond_eval_arg<std::false_type, a, Args...>
		{
			typedef a type;
		};

		template <template <class a, class b> class MetaFun, class a1, class a2>
		struct lambda_plExp<MetaFun<a1, a2> >
		{
			typedef struct
			{
				template <class A1 = void_, class A2 = void_>
				struct apply : MetaFun<
					typename cond_eval_arg<
						typename isPlaceholderExpression<a1>::type,
						a1,	A1, A2
					>::type,
					typename cond_eval_arg<
						typename isPlaceholderExpression<a2>::type,
						a2, A1, A2
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
					typename cond_eval_arg<
						typename isPlaceholderExpression<a1>::type,
						a1, A1, A2
					>::type,
					typename cond_eval_arg<
						typename isPlaceholderExpression<a2>::type,
						a2, A1, A2
					>::type,
					typename cond_eval_arg<
						typename isPlaceholderExpression<a3>::type,
						a3, A1, A2
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

	template <class Fun, class...Args>
	struct apply : lambda<Fun>::type::template apply<Args...> {};
}

#endif

