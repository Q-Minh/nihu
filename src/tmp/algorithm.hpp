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
 * \file algorithm.hpp
 * \ingroup tmp
 * \brief vectoralgorithms
 */

#ifndef ALGORITHM_HPP_INCLUDED
#define ALGORITHM_HPP_INCLUDED

#include <type_traits>

#include "integer.hpp"
#include "bool.hpp"
#include "lambda.hpp"
#include "sequence.hpp"
#include "operator.hpp"

namespace tmp
{
	namespace internal
	{
		/**
		 * \brief accumulate elements in a range using a user-specified metafunctor
		 * \tparam Beg begin iterator
		 * \tparam End end iterator
		 * \tparam Init initial value of accumulation
		 * \tparam Fun accumulating functor
		 */
		template <class Beg, class End, class Init, class Fun >
		struct accumulate_impl : accumulate_impl<
			typename next<Beg>::type,
			End,
			typename apply<Fun, Init, typename deref<Beg>::type>::type,
			Fun
		> {};

		/**
		 * \brief terminating case of accumulate_impl where the end iterator has been reached
		 */
		template <class End, class Init, class Fun>
		struct accumulate_impl<End, End, Init, Fun>
		{
			typedef Init type;
		};
	}

	/**
	 * \brief accumulate elements of a sequence using a binary metafunction
	 * \tparam Seq the sequence the elements of which are accumulated
	 * \tparam Init initial value of accumulation
	 * \tparam Fun accumulating functor, the default is plus
	 */
	template <class Seq, class Init, class Fun = tmp::plus<_1,_2> >
	struct accumulate : internal::accumulate_impl<
		typename begin<Seq>::type,
		typename end<Seq>::type,
		Init,
		Fun
	> {};

	/**
	 * \brief minimum of range
	 * \tparam Seq sequence
	 */
	template <class Seq>
	struct min : accumulate<
		Seq,
		typename deref<typename begin<Seq>::type>::type,
		if_<less<_1,_2>,_1,_2>
	> {};

	/**
	 * \brief maximum of range
	 * \tparam Seq sequence
	 */
	template <class Seq>
	struct max : accumulate<
		Seq,
		typename deref<typename begin<Seq>::type>::type,
		if_<less<_1,_2>,_2,_1>
	> {};


	namespace internal
	{
		class empty;

		template <class A, class B>
		struct inheriter
		{
			struct type : public A, public B {};
		};

		template <class B>
		struct inheriter<empty, B>
		{
			struct type : public B {};
		};
	}

	/**
	 * \brief combine a sequence of classes so that the result is inherited from each element
	 * \tparam Seq the sequence the elements of which are transformed
	 */
	template <class Seq>
	struct inherit : accumulate<
		Seq,
		internal::empty,
		internal::inheriter<_1,_2>
	> {};

	namespace internal
	{
		/**
		 * \brief transform elements in a range using a user-specified metafunctor and an inserter
		 * \tparam Beg begin iterator
		 * \tparam End end iterator
		 * \tparam Ins inserter used to fill output container
		 * \tparam Trans transformation functor
		 */
		template <class Beg, class End, class Ins, class Trans>
		struct transform_impl : transform_impl<
			typename next<Beg>::type,
			End,
			inserter<
				typename apply<
					typename Ins::operation,
					typename Ins::state,
					typename apply<
						Trans,
						typename deref<Beg>::type
					>::type
				>::type,
				typename Ins::operation
			>,
			Trans
		> {};

		/**
		 * \brief terminating case of transform_impl where the end iterator has been reached
		 */
		template <class End, class Ins, class Trans>
		struct transform_impl<End, End, Ins, Trans>
		{
			typedef typename Ins::state type;
		};

		/**
		 * \brief conditionally transform elements in a range
		 * \tparam Beg begin iterator
		 * \tparam End end iterator
		 * \tparam Ins inserter used to fill output container
		 * \tparam Cond condition evaluated for each element
		 * \tparam Trans transformation functor
		 */
		template <class Beg, class End, class Ins, class Cond, class Trans>
		struct transform_if_ptr_impl : transform_if_ptr_impl<
			typename next<Beg>::type,
			End,
			inserter<
				typename std::conditional<
					apply<Cond, Beg>::type::value,
					typename apply<
						typename Ins::operation,
						typename Ins::state,
						typename apply<Trans, Beg>::type
					>::type,
					typename Ins::state
				>::type,
				typename Ins::operation
			>,
			Cond,
			Trans
		> {};

		/**
		 * \brief terminating case of transform_impl where the end iterator has been reached
		 */
		template <class End, class Ins, class Cond, class Trans>
		struct transform_if_ptr_impl<End, End, Ins, Cond, Trans>
		{
			typedef typename Ins::state type;
		};
	}

	/**
	 * \brief transform elements in a sequence using a user-specified metafunctor and an inserter
	 * \tparam Seq the sequence the elements of which are transformed
	 * \tparam Ins inserter used to fill output container
	 * \tparam Trans transformation functor
	 */
	template <class Seq, class Ins, class Trans>
	struct transform : internal::transform_impl<
		typename begin<Seq>::type,
		typename end<Seq>::type,
		Ins,
		Trans
	> {};

	/**
	 * \brief conditionally transform elements in a sequence
	 * \tparam Seq the sequence the elements of which are transformed
	 * \tparam Ins inserter used to fill output container
	 * \tparam Cond condition
	 * \tparam Trans transformation functor
	 */
	template <class Seq, class Ins, class Cond, class Trans>
	struct transform_if_ptr : internal::transform_if_ptr_impl<
		typename begin<Seq>::type,
		typename end<Seq>::type,
		Ins,
		Cond,
		Trans
	> {};

	/**
	 * \brief copy elements from a container into an other
	 * \tparam Seq sequence
	 * \tparam Ins inserter used to fill output container
	 */
	template <class Seq, class Ins>
	struct copy : transform<Seq, Ins, _1> {};

	/**
	 * \brief copy elements from a container into an other
	 * \tparam Seq sequence
	 * \tparam Ins inserter used to fill output container
	 */
	template <class Seq, class Ins, class Cond>
	struct copy_if : transform_if_ptr<Seq, Ins, Cond, tmp::deref<_1> > {};

	/**
	 * \brief concatenate two sequences into a new sequence
	 * \tparam Seq1 first sequence
	 * \tparam Seq1 second sequence
	 */
	template <class Seq1, class Seq2>
	struct concatenate : transform<
		Seq2,
		inserter<Seq1, push_back<_1, _2> >,
		_1
	> {};

	/**
	 * \brief transform sequence of sequences into a flat sequence
	 * \tparam Seq sequence of sequences
	 */
	template <class Seq>
	struct serialise : accumulate<
		Seq,
		typename empty<Seq>::type,
		concatenate<_1, _2>
	> {};

	namespace internal
	{
		template <class Beg, class End, class Elem>
		struct find_impl : std::conditional<
			std::is_same<typename deref<Beg>::type, Elem>::type::value,
			Beg,
			typename find_impl<typename next<Beg>::type, End, Elem>::type
		> {};

		template <class End, class Elem>
		struct find_impl<End, End, Elem>
		{
			typedef End type;
		};
	}

	/**
	 * \brief Find an element in a sequence
	 * \details return an iterator for the first match. The end iterator is returned if the element is not found.
	 * \tparam Seq the sequence
	 * \tparam Elem the element
	 * \return iterator to the first occurrence of Elem in Seq or the end iterator if Elem is not found
	 */
	template <class Seq, class Elem>
	struct find : internal::find_impl<
		typename begin<Seq>::type,
		typename end<Seq>::type,
		Elem
	> {};

	/**
	 * \brief return true if the element is member of a sequence
	 * \details returns true_type if Elem is found in Seq
	 * \tparam Seq the sequence
	 * \tparam Elem the element
	 */
	template <class Seq, class Elem>
	struct is_member : not_<
		typename std::is_same<
			typename find<Seq, Elem>::type,
			typename end<Seq>::type
		>::type
	> {};

	/**
	 * \brief return a vector containing each element of Seq exactly once
	 * \tparam Seq a sequence
	 * \return a new sequence that contains each element of Seq exactly once
	 */
	template <class Seq>
	struct unique : accumulate<
		Seq,
		typename empty<Seq>::type,
		if_<is_member<_1,_2>, _1, push_back<_1,_2> >
	> {};

	namespace internal
	{
		/** \brief swap two elements if condition is true_type */
		template <class first, class second, class condition>
		struct swap_if;

		/** \brief specialisation of swap_if for the false case */
		template <class first, class second>
		struct swap_if<first, second, std::false_type>
		{
			typedef first first_t;
			typedef second second_t;
		};

		/** \brief specialisation of swap_if for the true case */
		template <class first, class second>
		struct swap_if<first, second, std::true_type>
		{
			typedef second first_t;
			typedef first second_t;
		};

		template <
			class Seq, class Compare,
			class siz = typename less<integer<int, 1>, typename size<Seq>::type>::type
		>
		struct bubble_cycle
		{
			// get first two elements and remove them from the sequence
			typedef typename deref<typename begin<Seq>::type>::type first_t;
			typedef typename deref<
				typename begin<typename pop_front<Seq>::type>::type
			>::type second_t;
			typedef typename pop_front<typename pop_front<Seq>::type>::type trunc;
			// swap the two elements if needed
			typedef typename apply<Compare, first_t, second_t>::type cond;
			typedef typename swap_if<first_t, second_t, cond>::first_t new_first;
			typedef typename swap_if<first_t, second_t, cond>::second_t new_second;
			// push the swapped second back
			typedef typename push_front<trunc, new_second>::type processed;
			// repeat for the remaining sequence and then push the first back
			typedef typename push_front<
				typename bubble_cycle<processed, Compare>::type,
				new_first
			>::type type;
		};

		template <class Seq, class Compare>
		struct bubble_cycle<Seq, Compare, std::false_type> : Seq {};
	}

	/** \brief sort a sequence by bubble sort
	 * \tparam Seq the sequence to sort
	 * \tparam Compare the comparison metfunctor
	 */
	template <
		class Seq,
		class Compare = less<_2, _1>,
		class cnt = typename size<Seq>::type
	>
	struct bubble_sort : bubble_sort<
		typename internal::bubble_cycle<Seq, Compare>::type,
		Compare,
		typename prev<cnt>::type
	> {};

	/** \brief terminating case of bubble_sort for short vectors */
	template <class Seq, class Compare>
	struct bubble_sort<Seq, Compare, integer<int, 0> > : Seq {};


	namespace internal
	{
		template <class value, unsigned N, class Seq>
		struct constant_sequence_impl : push_back<
			typename constant_sequence_impl<
				value, N-1, Seq
			>::type,
			value
		> {};

		template <class value, class Seq>
		struct constant_sequence_impl<value, 0, Seq>
		{
			typedef Seq type;
		};
	}

	/** \brief generate a constant sequence */
	template <class value, unsigned N, class Seq>
	struct constant_sequence : internal::constant_sequence_impl<
		value, N, typename empty<Seq>::type
	> {};
}

#endif // ALGORITHM_HPP_INCLUDED
