/**
 * \file sequence.hpp
 * \ingroup tmp
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 * \brief implementation of compile time sequences
 */
#ifndef SEQUENCE_HPP_INCLUDED
#define SEQUENCE_HPP_INCLUDED

namespace tmp
{
	namespace internal
	{
		template <class tag>
		struct empty_impl;

		template <class tag>
		struct size_impl;

		template <class tag>
		struct at_impl;

		template <class tag>
		struct begin_impl;

		template <class tag>
		struct end_impl;

		template <class tag>
		struct clear_impl;

		template <class tag>
		struct push_front_impl;

		template <class tag>
		struct push_back_impl;

		template <class tag>
		struct pop_front_impl;

		template <class tag>
		struct pop_back_impl;
	}

	/** \brief return empty sequence */
	template <class Seq>
	struct empty :
		internal::empty_impl<typename Seq::tag> {};

	/** \brief return size */
	template <class Seq>
	struct size :
		internal::size_impl<typename Seq::tag>::template apply<Seq> {};

	/** \brief return element at a given position */
	template <class Seq, class Pos>
	struct at :
		internal::at_impl<typename Seq::tag>::template apply<Seq, Pos> {};

	/** \brief return begin iterator of a sequence */
	template <class Seq>
	struct begin :
		internal::begin_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief return end iterator of a sequence */
	template <class Seq>
	struct end :
		internal::end_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief clear a sequence */
	template <class Seq>
	struct clear :
		internal::clear_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief push an element to the front */
	template <class Seq, class T>
	struct push_front :
		internal::push_front_impl<typename Seq::tag> :: template apply<Seq, T> {};

	/** \brief push an element to the back */
	template <class Seq, class T>
	struct push_back :
		internal::push_back_impl<typename Seq::tag> :: template apply<Seq, T> {};

	/** \brief pop an element from the front */
	template <class Seq>
	struct pop_front :
		internal::pop_front_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief pop an element from the back */
	template <class Seq>
	struct pop_back :
		internal::pop_back_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief metafunction to dereference an iterator */
	template <class Iter>
	struct deref;

	/** \brief a compile time inserter */
	template <class State, class Operation>
	struct inserter
	{
		/** \brief the actual state of the underlying sequence */
		typedef State state;
		/** \brief the operation performed by insertion */
		typedef Operation operation;
	};
}

#endif // SEQUENCE_HPP_INCLUDED

