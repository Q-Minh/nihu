#ifndef CALL_EACH_HPP
#define CALL_EACH_HPP

#include "sequence.hpp"
#include "lambda.hpp"
#include "algorithm.hpp"

template <class Begin, class End, class Transform, class Arg1, class Arg2>
struct call_each_str
{
	static void eval(Arg1 arg1, Arg2 arg2)
	{
		typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
		cur c;
		c(arg1, arg2);
		call_each_str<typename next<Begin>::type, End, Transform, Arg1, Arg2>::eval(arg1, arg2);
	}
};

template <class End, class Transform, class Arg1, class Arg2>
struct call_each_str<End, End, Transform, Arg1, Arg2>
{
	static void eval(Arg1, Arg2) { }
};


template <class Begin, class End, class Transform, class Arg1, class Arg2>
struct call_until_str
{
	static bool eval(Arg1 arg1, Arg2 arg2)
	{
		typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
		cur c;
		if (!c(arg1, arg2))
			return call_until_str<typename next<Begin>::type, End, Transform, Arg1, Arg2>::eval(arg1, arg2);
		return true;
	}
};

template <class End, class Transform, class Arg1, class Arg2>
struct call_until_str<End, End, Transform, Arg1, Arg2>
{
	static bool eval(Arg1, Arg2) { return false; }
};


template <class T>
struct is_reference : false_ {};

template <class T>
struct is_reference<T &> : true_ {};

template <class T>
struct add_ref
{
	typedef typename if_<
		typename is_reference<T>::type,
		T,
		T &
	>::type type;
};


template <class Seq, class Trans, class Arg1, class Arg2>
static void call_each(Arg1 &arg1, Arg2 &arg2)
{
	call_each_str<
		typename begin<Seq>::type,
		typename end<Seq>::type,
		Trans,
		typename add_ref<Arg1>::type,
		typename add_ref<Arg2>::type
	>::eval(arg1, arg2);
}


template <class Seq, class Trans, class Arg1, class Arg2>
static bool call_until(Arg1 &arg1, Arg2 &arg2)
{
	return call_until_str<
		typename begin<Seq>::type,
		typename end<Seq>::type,
		Trans,
		typename add_ref<Arg1>::type,
		typename add_ref<Arg2>::type
	>::eval(arg1, arg2);
}


#endif

