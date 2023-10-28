// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file lambda.hpp
 * \brief implementation of placeholders and lambda functions
 * \ingroup tmp
 */

#ifndef LAMBDA_HPP_INCLUDED
#define LAMBDA_HPP_INCLUDED

#include <type_traits>

namespace tmp
{
	/** \brief select N-th argument of a variadic template */
	template <unsigned N, class T, class ...Args>
	struct select_argument : select_argument<N-1, Args...> {};

	/** \brief terminating case of select_argument */
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

	/** \brief shorthand for selecting the 1st argument */
	typedef arg<1> _1;
	/** \brief shorthand for selecting the 2nd argument */
	typedef arg<2> _2;
	/** \brief shorthand for selecting the 3rd argument */
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

	// Forward declaration of the general case
	template <class T, class...Args>
	struct containsPlaceholderExpression;

	// General declaration of the inner impl class
	// Decision is made based on the first argument
	// (See specialisation below)
	template <class First, class...Args>
	struct containsPlaceholderExpressionImpl;

	// True branch
	// First parameter is a placeholder expression
	// Recursion tail for true_type
	template <class...Args>
	struct containsPlaceholderExpressionImpl<std::true_type, Args...> :
		std::true_type {};

	// Recursion tail for false_type
	template <>
	struct containsPlaceholderExpressionImpl<std::false_type> :
		std::false_type {};

	// General case, further recursion
	template <class...Args>
	struct containsPlaceholderExpressionImpl<std::false_type, Args...> :
		containsPlaceholderExpression<Args...> {};

	// General case
	template <class T, class...Args>
	struct containsPlaceholderExpression : containsPlaceholderExpressionImpl<
		typename isPlaceholderExpression<T>::type,
		Args...
	> {};

	// Specialisation for metafunction
	template <template <class...Args> class MF, class...Args>
	struct isPlaceholderExpression<MF<Args...> > :
		containsPlaceholderExpression<Args...> {};


	namespace internal
	{
		// The general case
		template <class Fun>
		struct lambda_plExp; // dummy case

		// Specialisation for placeholder
		template <unsigned N>
		struct lambda_plExp<arg<N> >
		{
			struct type
			{
				template <class...Args>
				struct apply : arg<N>::template apply<Args...> {};
			};
		};

		// Specialisation for metafunction for 1 parameter
		/** \todo variadic metafun would be nicer */
		template <template <class a> class MetaFun, class a1>
		struct lambda_plExp<MetaFun<a1> >
		{
			struct type
			{
				template <class...Args>
				struct apply
					: MetaFun<
						typename lambda_plExp<a1>::type::template apply<Args...>::type
						>
				{};
			};
		};

		// Conditional eveluation of arguments
		// First argument tells if a is a placeholder expression
		// General case
		template <class BOOL, class a, class...Args>
		struct cond_eval_arg;

		// Specialisation for true case (a is a plhexp)
		template <class a, class...Args>
		struct cond_eval_arg<std::true_type, a, Args...> :
			lambda_plExp<a>::type::template apply<Args...> {};

		// Specialisation for false case (a is not a plhexp)
		template <class a, class...Args>
		struct cond_eval_arg<std::false_type, a, Args...>
		{
			typedef a type;
		};

		// Two parameter metafun
		template <template <class a, class b> class MetaFun, class a1, class a2>
		struct lambda_plExp<MetaFun<a1, a2> >
		{
			struct type
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
					>::type
				> {};
			};
		};

		// Three parameter metafun (not used yet)
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
	struct lambda : std::conditional<
		isPlaceholderExpression<Exp>::type::value,
		typename internal::lambda_plExp<Exp>::type,
		Exp
	> {};

	/**
	 * \brief The apply metafunction shortcut for lambda evaluation
	 * \tparam Fun metafunction
	 * \tparam Args template arguments of Fun
	 */
	template <class Fun, class...Args>
	struct apply : lambda<Fun>::type::template apply<Args...> {};
}

#endif

