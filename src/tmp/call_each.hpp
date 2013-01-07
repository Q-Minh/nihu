#ifndef CALL_EACH_HPP
#define CALL_EACH_HPP

#include "sequence.hpp"
#include "lambda.hpp"

template <class Arg, class functor, class next>
struct inner_type 
{
	static void apply (Arg arg)
	{
		functor f;
		f(arg);

		next::apply(arg);
	}
};

template <class Arg>
struct inner_type_final
{
	static void apply (Arg) {}
};

template <class functor, class next>
struct inner_type<void, functor, next>
{
	static void apply (void)
	{
		functor f;
		f();

		next::apply();
	}
};

template <>
struct inner_type_final<void>
{
	static void apply (void) {}
};

template <class Begin, class End, class Fun, class Arg = void>
struct call_each
{
	typedef typename deref<Begin>::type obj;
	typedef typename lambda<Fun>::type::template apply<obj> functor;
	typedef typename call_each<typename next<Begin>::type, End, Fun, Arg>::type next;
	
	struct type : inner_type<Arg, functor, next> {};
};

template <class End, class Fun, class Arg>
struct call_each<End, End, Fun, Arg>
{
	struct type : inner_type_final<Arg> {};
};

#endif

