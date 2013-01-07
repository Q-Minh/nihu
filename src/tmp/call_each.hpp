#ifndef CALL_EACH_HPP
#define CALL_EACH_HPP

#include "sequence.hpp"
#include "lambda.hpp"

template <class functor, class next, class Arg1, class Arg2>
struct inner_type 
{
	static void apply (Arg1 arg1, Arg2 arg2)
	{
		functor f;
		f(arg1, arg2);

		next::apply(arg1, arg2);
	}
};

template <class functor, class next, class Arg1>
struct inner_type <functor, next, Arg1, void>
{
	static void apply (Arg1 arg1)
	{
		functor f;
		f(arg1);

		next::apply(arg1);
	}
};

template <class functor, class next>
struct inner_type<functor, next, void, void>
{
	static void apply (void)
	{
		functor f;
		f();

		next::apply();
	}
};

template <class Arg1, class Arg2>
struct inner_type_final
{
	static void apply (Arg1, Arg2) {}
};

template <class Arg1>
struct inner_type_final<Arg1, void>
{
	static void apply (Arg1) {}
};

template <>
struct inner_type_final<void, void>
{
	static void apply (void) {}
};

template <class Begin, class End, class Fun, class Arg1 = void, class Arg2 = void>
struct call_each
{
	typedef typename deref<Begin>::type obj;
	typedef typename lambda<Fun>::type::template apply<obj> functor;
	typedef typename call_each<typename next<Begin>::type, End, Fun, Arg1, Arg2>::type next;
	
	typedef inner_type<functor, next, Arg1, Arg2> type;
};

template <class End, class Fun, class Arg1, class Arg2>
struct call_each<End, End, Fun, Arg1, Arg2>
{
	typedef inner_type_final<Arg1, Arg2> type;
};


#endif
