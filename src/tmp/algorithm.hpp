/**
 * \file algorithm.hpp
 * \brief vectoralgorithms
 */

#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "bool.hpp"
#include "general.hpp"
#include "lambda.hpp"
#include "sequence.hpp"
#include "vector.hpp"
#include "operator.hpp"

namespace tmp
{
	namespace internal
	{
		/**
		 * \brief accumulate_impl elements in a range using a user-specified metafunctor
		 * \tparam Beg begin iterator
		 * \tparam End end iterator
		 * \tparam Init initial value of accumulation
		 * \tparam Fun accumulating functor, the default is plus
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


	class empty;

	namespace internal
	{
		template <class A, class B>
		struct inheriter
		{
			struct type : public A, B {};
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
	struct inherit : accumulate<Seq, empty, internal::inheriter<_1,_2> > {};
	
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
	 * \brief copy elements from a container into an other
	 * \tparam Seq sequence
	 * \tparam Ins inserter used to fill output container
	 */
	template <class Seq, class Ins>
	struct copy : transform<Seq, Ins, _1> {};

	/**
	 * \brief concatenate two sequences into a tiny sequence
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
		vector<>,
		concatenate<_1, _2>
	> {};

	namespace internal
	{
		template <class Beg, class End, class Elem>
		struct find_impl : if_<
			typename is_same<typename deref<Beg>::type, Elem>::type,
			Beg,
			typename find_impl<typename next<Beg>::type, End, Elem>::type
		> {};

		template <class End, class Elem>
		struct find_impl<End, End, Elem>
		{
			typedef End type;
		};
	}

	template <class Seq, class Elem>
	struct find : internal::find_impl<
		typename begin<Seq>::type,
		typename end<Seq>::type,
		Elem
	> {};

	template <class Seq, class Elem>
	struct is_member : not_<
		typename is_same<
			typename find<Seq, Elem>::type,
			typename end<Seq>::type
		>::type
	> {};

	template <class Seq>
	struct unique : accumulate<
		Seq,
		vector<>,
		if_<is_member<_1,_2>, _1, push_back<_1,_2> >
	> {};
}

#endif

