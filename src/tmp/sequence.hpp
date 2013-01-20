/**
 * \file sequence.hpp
 * \brief implementation of a vector
 */
#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include "integer.hpp"

namespace tmp
{
	namespace internal
	{
		/**
		 * \brief metafunctor to implement size operation
		 */
		template <class tag>
		struct size_impl;

		/**
		 * \brief metafunctor to implement at operation
		 */
		template <class tag>
		struct at_impl;

		/**
		 * \brief metafunctor to implement begin operation
		 */
		template <class tag>
		struct begin_impl;

		/**
		 * \brief metafunctor to implement end operation
		 */
		template <class tag>
		struct end_impl;

		/**
		 * \brief metafunctor to implement clear operation
		 */
		template <class tag>
		struct clear_impl;

		/**
		 * \brief metafunctor to implement push_front operation
		 */
		template <class tag>
		struct push_front_impl;

		/**
		 * \brief metafunctor to implement push_back operation
		 */
		template <class tag>
		struct push_back_impl;

		/**
		 * \brief metafunctor to implement pop_back operation
		 */
		template <class tag>
		struct pop_back_impl;
	}

	/**
	 * \brief metafunction returning size (uses size_impl metafunctor)
	 */
	template <class Seq>
	struct size : internal::size_impl<typename Seq::tag>::template apply<Seq> {};

	/**
	 * \brief metafunction returning element at a given position (uses at_impl metafunctor)
	 */
	template <class Seq, class Pos>
	struct at : internal::at_impl<typename Seq::tag>::template apply<Seq, Pos> {};

	/**
	 * \brief metafunction returning begin iterator of a sequence (uses begin_impl metafunctor)
	 */
	template <class Seq>
	struct begin : internal::begin_impl<typename Seq::tag> :: template apply<Seq> {};

	/**
	 * \brief metafunction returning end iterator of a sequence (uses end_impl metafunctor)
	 */
	template <class Seq>
	struct end : internal::end_impl<typename Seq::tag> :: template apply<Seq> {};

	/**
	 * \brief metafunction clearing a sequence (uses clear_impl metafunctor)
	 */
	template <class Seq>
	struct clear : internal::clear_impl<typename Seq::tag> :: template apply<Seq> {};

	/**
	 * \brief metafunction pushing an element to the front (uses push_front_impl metafunctor)
	 */
	template <class Seq, class T>
	struct push_front : internal::push_front_impl<typename Seq::tag> :: template apply<Seq, T> {};

	/**
	 * \brief metafunction pushing an element to the back (uses push_back_impl metafunctor)
	 */
	template <class Seq, class T>
	struct push_back : internal::push_back_impl<typename Seq::tag> :: template apply<Seq, T> {};

	/**
	 * \brief metafunction popping an element from the back (uses pop_back_impl metafunctor)
	 */
	template <class Seq>
	struct pop_back : internal::pop_back_impl<typename Seq::tag> :: template apply<Seq> {};

	template <class Iter>
	struct deref;

	template <class State, class Operation>
	struct inserter
	{
		typedef State state;
		typedef Operation operation;
	};
}

#endif

