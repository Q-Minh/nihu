/**
 * \file control.hpp
 * \author Peter fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 * \brief Implementation of code generating control structures
 */
#ifndef CALL_EACH_HPP
#define CALL_EACH_HPP

#include <type_traits>

#include "sequence.hpp"
#include "lambda.hpp"
#include "algorithm.hpp"

/**
 * \brief template metaprogramming functions
 */
namespace tmp
{
	namespace internal
	{
		template <class Begin, class End, class Transform, class Arg1 = void, class Arg2 = void>
		struct call_each_impl
		{
			static void eval(Arg1 arg1, Arg2 arg2)
			{
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				c(arg1, arg2);
				call_each_impl<typename next<Begin>::type, End, Transform, Arg1, Arg2>::eval(arg1, arg2);
			}
		};

		template <class End, class Transform, class Arg1, class Arg2>
		struct call_each_impl<End, End, Transform, Arg1, Arg2>
		{
			static void eval(Arg1, Arg2) { }
		};


		template <class Begin, class End, class Transform, class Arg1>
		struct call_each_impl<Begin, End, Transform, Arg1, void>
		{
			static void eval(Arg1 arg1)
			{
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				c(arg1);
				call_each_impl<typename next<Begin>::type, End, Transform, Arg1>::eval(arg1);
			}
		};

		template <class End, class Transform, class Arg1>
		struct call_each_impl<End, End, Transform, Arg1, void>
		{
			static void eval(Arg1) { }
		};


		template <class Begin, class End, class Transform>
		struct call_each_impl<Begin, End, Transform, void, void>
		{
			static void eval(void)
			{
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				c();
				call_each_impl<typename next<Begin>::type, End, Transform>::eval();
			}
		};

		template <class End, class Transform>
		struct call_each_impl<End, End, Transform, void, void>
		{
			static void eval(void) { }
		};
	}

	/**
	 * \brief call function object with different types
	 * \tparam Seq a type sequence
	 * \tparam Trans a function object
	 */
	template <class Seq, class Transform>
	static void call_each(void)
	{
		internal::call_each_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform
		>::eval();
	}

	/** \brief one argument version of call_each */
	template <class Seq, class Transform, class Arg1>
	static void call_each(Arg1 arg1)
	{
		internal::call_each_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform,
			typename std::add_lvalue_reference<Arg1>::type
		>::eval(arg1);
	}


	/** \brief two argument version of call_each */
	template <class Seq, class Transform, class Arg1, class Arg2>
	static void call_each(Arg1 arg1, Arg2 arg2)
	{
		internal::call_each_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform,
			typename std::add_lvalue_reference<Arg1>::type,
			typename std::add_lvalue_reference<Arg2>::type
		>::eval(arg1, arg2);
	}

	namespace internal
	{
		template <class Begin, class End, class SeqIn, class Transform, class Arg1 = void, class Arg2 = void>
		struct d_call_each_impl
		{
			typedef typename lambda<Transform>::type::template apply<
				typename deref<Begin>::type, _1
			> partially_evaluated_transform;
			static void eval(Arg1 arg1, Arg2 arg2)
			{
				call_each<
					SeqIn,
					partially_evaluated_transform,
					Arg1,
					Arg2
				>(arg1, arg2);
				d_call_each_impl<typename next<Begin>::type, End, SeqIn, Transform, Arg1, Arg2>::eval(arg1, arg2);
			}
		};


		template <class End, class SeqIn, class Transform, class Arg1, class Arg2>
		struct d_call_each_impl<End, End, SeqIn, Transform, Arg1, Arg2>
		{
			static void eval(Arg1 arg1, Arg2 arg2) {}
		};
	}


	template <class SeqOut, class SeqIn, class Trans, class Arg1, class Arg2>
	static void d_call_each(Arg1 a1, Arg2 a2)
	{
		internal::d_call_each_impl<
					typename begin<SeqOut>::type,
					typename end<SeqOut>::type,
					SeqIn,
					Trans,
					typename std::add_lvalue_reference<Arg1>::type,
					typename std::add_lvalue_reference<Arg2>::type
				>::eval(a1, a2);
	}


	namespace internal
	{
		template <class Begin, class End, class Transform, class Arg1 = void, class Arg2 = void>
		struct call_until_impl
		{
			static bool eval(Arg1 arg1, Arg2 arg2)
			{
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				if (!c(arg1, arg2))
					return call_until_impl<typename next<Begin>::type, End, Transform, Arg1, Arg2>::eval(arg1, arg2);
				return true;
			}
		};

		template <class End, class Transform, class Arg1, class Arg2>
		struct call_until_impl<End, End, Transform, Arg1, Arg2>
		{
			static bool eval(Arg1, Arg2) { return false; }
		};


		template <class Begin, class End, class Transform, class Arg1>
		struct call_until_impl<Begin, End, Transform, Arg1, void>
		{
			static bool eval(Arg1 arg1)
			{
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				if (!c(arg1))
					return call_until_impl<typename next<Begin>::type, End, Transform, Arg1>::eval(arg1);
				return true;
			}
		};

		template <class End, class Transform, class Arg1>
		struct call_until_impl<End, End, Transform, Arg1, void>
		{
			static bool eval(Arg1) { return false; }
		};


		template <class Begin, class End, class Transform>
		struct call_until_impl<Begin, End, Transform, void, void>
		{
			static bool eval(void)
			{
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				if (!c())
					return call_until_impl<typename next<Begin>::type, End, Transform>::eval();
				return true;
			}
		};

		template <class End, class Transform>
		struct call_until_impl<End, End, Transform, void, void>
		{
			static bool eval(void) { return false; }
		};
	}

	/**
	 * \brief call function object with different types until it returns true
	 * \tparam Seq a type sequence
	 * \tparam Trans a function object
	 */
	template <class Seq, class Transform>
	static bool call_until(void)
	{
		return internal::call_until_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform
		>::eval();
	}


	/** \brief one argument version of call_until */
	template <class Seq, class Transform, class Arg1>
	static bool call_until(Arg1 &arg1)
	{
		return internal::call_until_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform,
			typename std::add_lvalue_reference<Arg1>::type
		>::eval(arg1);
	}


	/** \brief two argument version of call_until */
	template <class Seq, class Transform, class Arg1, class Arg2>
	static bool call_until(Arg1 &arg1, Arg2 &arg2)
	{
		return internal::call_until_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform,
			typename std::add_lvalue_reference<Arg1>::type,
			typename std::add_lvalue_reference<Arg2>::type
		>::eval(arg1, arg2);
	}
}

#endif

