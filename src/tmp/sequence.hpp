/**
 * \file sequence.hpp
 * \ingroup tmp
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu 
 * \brief implementation of compile time sequences
 */
#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include "integer.hpp"

namespace tmp
{
	namespace internal
	{
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

	/** \brief metafunction returning size */
	template <class Seq>
	struct size :
		internal::size_impl<typename Seq::tag>::template apply<Seq> {};

	/** \brief metafunction returning element at a given position */
	template <class Seq, class Pos>
	struct at :
		internal::at_impl<typename Seq::tag>::template apply<Seq, Pos> {};

	/** \brief metafunction returning begin iterator of a sequence */
	template <class Seq>
	struct begin :
		internal::begin_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief metafunction returning end iterator of a sequence */
	template <class Seq>
	struct end :
		internal::end_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief metafunction clearing a sequence */
	template <class Seq>
	struct clear :
		internal::clear_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief metafunction pushing an element to the front */
	template <class Seq, class T>
	struct push_front :
		internal::push_front_impl<typename Seq::tag> :: template apply<Seq, T> {};

	/** \brief metafunction pushing an element to the back */
	template <class Seq, class T>
	struct push_back :
		internal::push_back_impl<typename Seq::tag> :: template apply<Seq, T> {};

	/** \brief metafunction popping an element from the front */
	template <class Seq>
	struct pop_front :
		internal::pop_front_impl<typename Seq::tag> :: template apply<Seq> {};

	/** \brief metafunction popping an element from the back */
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
		typedef State state;
		typedef Operation operation;
	};
}

#endif

