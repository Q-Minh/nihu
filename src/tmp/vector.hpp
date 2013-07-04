/**
 * \file vector.hpp
 * \ingroup tmp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of a compile time vector ::tmp::vector
 */
#ifndef VECTOR_HPP_INCLUDED
#define VECTOR_HPP_INCLUDED

#include "sequence.hpp"
#include <tuple>

namespace tmp
{
	/** \brief the vector tag used for the selection of algorithms */
	struct vector_tag;

	/** \brief a compile time vector with an arbitrary number of arguments
	* \details the implementation is now based on std::tuple
	*/
	template <class...Args>
	struct vector
	{
		/** \brief self returning */
		typedef vector type;
		/** \brief the tag for algorithm selection */
		typedef vector_tag tag;

		/** \brief the underlying std::tuple container */
		typedef std::tuple<Args...> impl;
	};

	namespace internal
	{
		template <>
		struct at_impl<vector_tag>
		{
			template <class Seq, class Pos>
			struct apply : std::tuple_element<Pos::value, typename Seq::impl> {};
		};
	}

	/** \brief the vector iterator type used by tmp::begin, tmp::end and tmp::deref */
	template <class Seq, class Pos>
	struct vector_iterator;

	/** \brief increment a vector iterator */
	template <class Seq, class Pos>
	struct next<vector_iterator<Seq, Pos> >
	{
		typedef vector_iterator<Seq, typename next<Pos>::type> type;
	};

	/** \brief decrement a vector iterator */
	template <class Seq, class Pos>
	struct prev<vector_iterator<Seq, Pos> >
	{
		typedef vector_iterator<Seq, typename prev<Pos>::type> type;
	};

	namespace internal
	{
		template <>
		struct size_impl<vector_tag>
		{
			template <class Seq>
			struct apply : int_<std::tuple_size<typename Seq::impl>::value> {};
		};
	}

	/** \brief dereference a vector iterator */
	template <class Seq, class Pos>
	struct deref<vector_iterator<Seq, Pos> > : at<Seq, Pos> {};

	namespace internal
	{
		template <>
		struct begin_impl<vector_tag>
		{
			template <class Vect>
			struct apply
			{
				typedef vector_iterator<Vect, int_<0> > type;
			};
		};

		template <>
		struct end_impl<vector_tag>
		{
			template <class Vect>
			struct apply
			{
				typedef vector_iterator<
					Vect,
					typename size_impl<vector_tag>::template apply<Vect>::type
				> type;
			};
		};

		template <>
		struct clear_impl<vector_tag>
		{
			template <class Seq>
			struct apply : vector<> {};
		};

		/** \brief push a new type to the front of a tmp::vector */
		template <class Vect, class T>
		struct vector_push_front;

		template <class...Args, class T>
		struct vector_push_front<vector<Args...>, T> : vector<T, Args...> {};

		template <>
		struct push_front_impl<vector_tag>
		{
			template <class Vect, class T>
			struct apply : vector_push_front<Vect, T> {};
		};

		/** \brief push a new type to the back of a tmp::vector */
		template <class Vect, class T>
		struct vector_push_back;

		template <class...Args, class T>
		struct vector_push_back<vector<Args...>, T> : vector<Args..., T> {};

		template <>
		struct push_back_impl<vector_tag>
		{
			template <class Seq, class T>
			struct apply : vector_push_back<Seq, T> {};
		};


		/** \brief pop a type from the front of a tmp::vector */
		template <class Vect>
		struct vector_pop_front;

		template <class...Args, class T>
		struct vector_pop_front<vector<T, Args...> > : vector<Args...> {};

		template <>
		struct pop_front_impl<vector_tag>
		{
			template <class Seq>
			struct apply : vector_pop_front<Seq> {};
		};


/*
		template <class Vect>
		struct vector_pop_back;

		template <class...Args, class T>
		struct vector_pop_back<vector<Args..., T> > : vector<Args...> {};

		template <>
		struct pop_back_impl<vector_tag>
		{
			template <class Seq>
			struct apply : vector_pop_back<Seq> {};
		};
*/
	}
}

#endif

