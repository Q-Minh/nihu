/**
 * \file vector.hpp
 * \brief implementation of tmp::vector
 */
#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "sequence.hpp"

namespace tmp
{
	namespace internal
	{
		struct none;
		struct vector_tag;
	}

	template <class A0 = internal::none, class A1 = internal::none, class A2 = internal::none, class A3 = internal::none>
	struct vector
	{
		typedef vector type; 	/* self-returning */

		typedef internal::vector_tag tag;	/* tagged */

		typedef A0 arg0;
		typedef A1 arg1;
		typedef A2 arg2;
		typedef A3 arg3;
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

		template <class Seq, class Pos>
		struct deref<vector_iterator<Seq, Pos> > : at<Seq, Pos> {};

		template <>
		struct begin_impl<vector_tag>
		{
			template <class vector>
			struct apply
			{ /* metafunction forwarding impossible because vector_iterator is incomplete type */
				typedef vector_iterator<vector, int_<0> > type;
			};
		};

		template <>
		struct end_impl<vector_tag>
		{
			template <class vector>
			struct apply
			{
				typedef vector_iterator<
					vector,
					typename vector_size<typename vector::arg0, typename vector::arg1, typename vector::arg2, typename vector::arg3>::type
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
			template <class vector, class T>
			struct apply : vector<T, typename vector::arg0, typename vector::arg1, typename vector::arg2> {};
		};

		template <class vector, class T, int N>
		struct vector_push_back;

		template <class vector, class T>
		struct vector_push_back<vector, T, 0> : vector<T> {};

		template <class vector, class T>
		struct vector_push_back<vector, T, 1> : vector<typename vector::arg0, T> {};

		template <class vector, class T>
		struct vector_push_back<vector, T, 2> : vector<typename vector::arg0, typename vector::arg1, T> {};

		template <class vector, class T>
		struct vector_push_back<vector, T, 3> : vector<typename vector::arg0, typename vector::arg1, typename vector::arg2, T> {};

		template <>
		struct push_back_impl<vector_tag>
		{
			template <class Seq, class T>
			struct apply : vector_push_back<Seq, T, size<Seq>::value > {};
		};

		template <class vector, int N>
		struct vector_pop_back;

		template <class vector>
		struct vector_pop_back<vector, 1> : vector<> {};

		template <class vector>
		struct vector_pop_back<vector, 2> : vector<typename vector::arg0> {};

		template <class vector>
		struct vector_pop_back<vector, 3> : vector<typename vector::arg0, typename vector::arg1> {};

		template <class vector>
		struct vector_pop_back<vector, 4> : vector<typename vector::arg0, typename vector::arg1, typename vector::arg2> {};

		template <>
		struct pop_back_impl<vector_tag>
		{
			template <class Seq>
			struct apply : vector_pop_back<Seq, size<Seq>::value > {};
		};
	}
}

#endif

