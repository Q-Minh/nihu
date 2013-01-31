/**
 * \file vector.hpp
 * \brief implementation of tmp::vector
 */
#ifndef VECTOR_HPP_INCLUDED
#define VECTOR_HPP_INCLUDED

#include "sequence.hpp"

namespace tmp
{
	struct none;
	struct vector_tag;

	/**
	 * \brief A random access vector with up to 4 elements
	 * \tparam A0 0-th element
	 * \tparam A1 1-st element
	 * \tparam A2 2-nd element
	 * \tparam A3 3-rd element
	 */
	template <class A0 = none, class A1 = none, class A2 = none, class A3 = none>
	struct vector
	{
		typedef vector type; 	/**< \brief self-returning */

		typedef vector_tag tag;	/**< \brief tagged */

		typedef A0 arg0;	/**< \brief 0-th element */
		typedef A1 arg1;	/**< \brief 1-th element */
		typedef A2 arg2;	/**< \brief 2-th element */
		typedef A3 arg3;	/**< \brief 3-th element */
	};

	namespace internal
	{
		template <class Arg, int Pos>
		struct vector_at;

		template <class A0, class A1, class A2, class A3>
		struct vector_at<vector<A0, A1, A2, A3>, 0> { typedef A0 type; };

		template <class A0, class A1, class A2, class A3>
		struct vector_at<vector<A0, A1, A2, A3>, 1> { typedef A1 type; };

		template <class A0, class A1, class A2, class A3>
		struct vector_at<vector<A0, A1, A2, A3>, 2> { typedef A2 type; };

		template <class A0, class A1, class A2, class A3>
		struct vector_at<vector<A0, A1, A2, A3>, 3> { typedef A3 type; };

		template <>
		struct at_impl<vector_tag>
		{
			template <class Seq, class Pos>
			struct apply : vector_at<Seq, Pos::value> {};
		};
	}

	template <class Seq, class Pos>
	struct vector_iterator;

	template <class Seq, class Pos>
	struct next<vector_iterator<Seq, Pos> >
	{ /* metafunction forwarding impossible because vector_iterator is incomplete type */
		typedef vector_iterator<Seq, typename next<Pos>::type> type;
	};

	template <class Seq, class Pos>
	struct prev<vector_iterator<Seq, Pos> >
	{ /* metafunction forwarding impossible because vector_iterator is incomplete type */
		typedef vector_iterator<Seq, typename prev<Pos>::type> type;
	};

	namespace internal
	{
		template <class T0, class T1, class T2, class T3>
		struct vector_size : int_<4> {};

		template <class T0, class T1, class T2>
		struct vector_size<T0, T1, T2, none> : int_<3> {};

		template <class T0, class T1>
		struct vector_size<T0, T1, none, none> : int_<2> {};

		template <class T0>
		struct vector_size<T0, none, none, none> : int_<1> {};

		template <>
		struct vector_size<none, none, none, none> : int_<0> {};

		template <>
		struct size_impl<vector_tag>
		{
			template <class Seq>
			struct apply : vector_size<typename Seq::arg0, typename Seq::arg1, typename Seq::arg2, typename Seq::arg3> {};
		};
	}

	template <class Seq, class Pos>
	struct deref<vector_iterator<Seq, Pos> > : at<Seq, Pos> {};

	namespace internal
	{
		template <>
		struct begin_impl<vector_tag>
		{
			template <class Vect>
			struct apply
			{ /* metafunction forwarding impossible because vector_iterator is incomplete type */
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
					typename vector_size<typename Vect::arg0, typename Vect::arg1, typename Vect::arg2, typename Vect::arg3>::type
				> type;
			};
		};

		template <>
		struct clear_impl<vector_tag>
		{
			template <class Seq>
			struct apply : vector<> {};
		};

		template <>
		struct push_front_impl<vector_tag>
		{
			template <class Vect, class T>
			struct apply : vector<T, typename Vect::arg0, typename Vect::arg1, typename Vect::arg2> {};
		};

		template <class Vect, class T, int N>
		struct vector_push_back;

		template <class Vect, class T>
		struct vector_push_back<Vect, T, 0> : vector<T> {};

		template <class Vect, class T>
		struct vector_push_back<Vect, T, 1> : vector<typename Vect::arg0, T> {};

		template <class Vect, class T>
		struct vector_push_back<Vect, T, 2> : vector<typename Vect::arg0, typename Vect::arg1, T> {};

		template <class Vect, class T>
		struct vector_push_back<Vect, T, 3> : vector<typename Vect::arg0, typename Vect::arg1, typename Vect::arg2, T> {};

		template <>
		struct push_back_impl<vector_tag>
		{
			template <class Seq, class T>
			struct apply : vector_push_back<Seq, T, size<Seq>::value > {};
		};

		template <>
		struct pop_front_impl<vector_tag>
		{
			template <class Vect>
			struct apply : vector<typename Vect::arg0, typename Vect::arg1, typename Vect::arg2> {};
		};

		template <class Vect, int N>
		struct vector_pop_back;

		template <class Vect>
		struct vector_pop_back<Vect, 1> : vector<> {};

		template <class Vect>
		struct vector_pop_back<Vect, 2> : vector<typename Vect::arg0> {};

		template <class Vect>
		struct vector_pop_back<Vect, 3> : vector<typename Vect::arg0, typename Vect::arg1> {};

		template <class Vect>
		struct vector_pop_back<Vect, 4> : vector<typename Vect::arg0, typename Vect::arg1, typename Vect::arg2> {};

		template <>
		struct pop_back_impl<vector_tag>
		{
			template <class Vect>
			struct apply : vector_pop_back<Vect, size<Vect>::value > {};
		};
	}
}

#endif

