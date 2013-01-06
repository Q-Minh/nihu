#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

/**
 * \file sequence.hpp
 * \brief implementation of a vector
 */
#include "integer.hpp"

/**
 * \brief metafunctor to implement size operation
 */
template <class tag>
struct size_impl;

/**
 * \brief metafunction returning size (uses size_impl metafunctor)
 */
template <class Seq>
struct size : size_impl<typename Seq::tag>::template apply<Seq> {};

/**
 * \brief metafunctor to implement at operation
 */
template <class tag>
struct at_impl;

/**
 * \brief metafunction returning element at a given position (uses at_impl metafunctor)
 */
template <class Seq, class Pos>
struct at : at_impl<typename Seq::tag>::template apply<Seq, Pos> {};

/**
 * \brief metafunctor to implement begin operation
 */
template <class tag>
struct begin_impl;

/**
 * \brief metafunction returning begin iterator of a sequence (uses begin_impl metafunctor)
 */
template <class Seq>
struct begin : begin_impl<typename Seq::tag> :: template apply<Seq> {};

/**
 * \brief metafunctor to implement end operation
 */
template <class tag>
struct end_impl;

/**
 * \brief metafunction returning end iterator of a sequence (uses end_impl metafunctor)
 */
template <class Seq>
struct end : end_impl<typename Seq::tag> :: template apply<Seq> {};

/**
 * \brief metafunctor to implement clear operation
 */
template <class tag>
struct clear_impl;

/**
 * \brief metafunction clearing a sequence (uses clear_impl metafunctor)
 */
template <class Seq>
struct clear : clear_impl<typename Seq::tag> :: template apply<Seq> {};

/**
 * \brief metafunctor to implement push_front operation
 */
template <class tag>
struct push_front_impl;

/**
 * \brief metafunction pushing an element to the front (uses push_front_impl metafunctor)
 */
template <class Seq, class T>
struct push_front : push_front_impl<typename Seq::tag> :: template apply<Seq, T> {};

/**
 * \brief metafunctor to implement push_back operation
 */
template <class tag>
struct push_back_impl;

/**
 * \brief metafunction pushing an element to the back (uses push_back_impl metafunctor)
 */
template <class Seq, class T>
struct push_back : push_back_impl<typename Seq::tag> :: template apply<Seq, T> {};


template <class Iter>
struct deref;

struct none;

struct tiny_tag;

template <class A0 = none, class A1 = none, class A2 = none>
struct tiny
{
	typedef tiny type; 	/* self-returning */

	typedef tiny_tag tag;	/* tagged */

	typedef A0 arg0;
	typedef A1 arg1;
	typedef A2 arg2;
};


template <class Arg, int Pos>
struct tiny_at;

template <class A0, class A1, class A2>
struct tiny_at<tiny<A0, A1, A2>, 0> { typedef A0 type; };

template <class A0, class A1, class A2>
struct tiny_at<tiny<A0, A1, A2>, 1> { typedef A1 type; };

template <class A0, class A1, class A2>
struct tiny_at<tiny<A0, A1, A2>, 2> { typedef A2 type; };

template <>
struct at_impl<tiny_tag>
{
	template <class Seq, class Pos>
	struct apply : tiny_at<Seq, Pos::value> {};
};


template <class Seq, class Pos>
struct tiny_iterator;

template <class Seq, class Pos>
struct next<tiny_iterator<Seq, Pos> >
{ /* metafunction forwarding impossible because tiny_iterator is incomplete type */
	typedef tiny_iterator<Seq, typename next<Pos>::type> type;
};

template <class Seq, class Pos>
struct prev<tiny_iterator<Seq, Pos> >
{ /* metafunction forwarding impossible because tiny_iterator is incomplete type */
	typedef tiny_iterator<Seq, typename prev<Pos>::type> type;
};

template <class T0, class T1, class T2>
struct tiny_size : int_<3> {};

template <class T0, class T1>
struct tiny_size<T0, T1, none> : int_<2> {};

template <class T0>
struct tiny_size<T0, none, none> : int_<1> {};

template <>
struct tiny_size<none, none, none> : int_<0> {};

template <>
struct size_impl<tiny_tag>
{
	template <class Seq>
	struct apply : tiny_size<typename Seq::arg0, typename Seq::arg1, typename Seq::arg2> {};
};

template <class Seq, class Pos>
struct deref<tiny_iterator<Seq, Pos> > : at<Seq, Pos> {};

template <>
struct begin_impl<tiny_tag>
{
	template <class Tiny>
	struct apply
	{ /* metafunction forwarding impossible because tiny_iterator is incomplete type */
		typedef tiny_iterator<Tiny, int_<0> > type;
	};
};

template <>
struct end_impl<tiny_tag>
{
	template <class Tiny>
	struct apply
	{
		typedef tiny_iterator<
			Tiny,
			typename tiny_size<typename Tiny::arg0, typename Tiny::arg1, typename Tiny::arg2>::type
		> type;
	};
};

template <>
struct clear_impl<tiny_tag>
{
	template <class Seq>
	struct apply : tiny<> {};
};

template <>
struct push_front_impl<tiny_tag>
{
	template <class Tiny, class T>
	struct apply : tiny<T, typename Tiny::arg0, typename Tiny::arg1> {};
};

template <class Tiny, class T, int N>
struct tiny_push_back;

template <class Tiny, class T>
struct tiny_push_back<Tiny, T, 0> : tiny<T> {};

template <class Tiny, class T>
struct tiny_push_back<Tiny, T, 1> : tiny<typename Tiny::arg0, T> {};

template <class Tiny, class T>
struct tiny_push_back<Tiny, T, 2> : tiny<typename Tiny::arg0, typename Tiny::arg1, T> {};

template <>
struct push_back_impl<tiny_tag>
{
	template <class Seq, class T>
	struct apply : tiny_push_back<Seq, T, size<Seq>::value > {};
};

#endif

