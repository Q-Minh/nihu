/**
 * \file stack.hpp
 * \brief implementation of tmp::stack
 */
#ifndef STACK_HPP_INCLUDED
#define STACK_HPP_INCLUDED

#include "integer.hpp"
#include "sequence.hpp"

namespace tmp
{
	struct none;
	struct stack_tag;

	namespace internal
	{
		template <class Head = none, class Tail = none>
		struct stack_impl
		{
			typedef stack_impl type;

			typedef stack_tag tag;

			typedef Head head;
			typedef Tail tail;
		};
	}
	
	/** \brief an empty stack */
	typedef internal::stack_impl<none, none> empty_stack;

	template <class Seq>
	struct stack_iterator;

	template <class Seq>
	struct next<stack_iterator<Seq> >
	{
		typedef stack_iterator<typename Seq::tail> type;
	};

	namespace internal
	{
		template <class Stack>
		struct stack_size : int_<stack_size<typename Stack::tail>::type::value + 1> {};

		template <>
		struct stack_size<empty_stack> : int_<0> {};

		template <>
		struct size_impl<stack_tag>
		{
			template <class Seq>
			struct apply : stack_size<Seq> {};
		};
	}

	/** \brief specialisation of metafunction ::deref for the stack iterator */
	template <class Seq>
	struct deref<stack_iterator<Seq> >
	{
		typedef typename Seq::head type;
	};

	namespace internal
	{
		template <>
		struct begin_impl<stack_tag>
		{
			template <class Stack>
			struct apply
			{
				typedef stack_iterator<Stack> type;
			};
		};

		template <>
		struct end_impl<stack_tag>
		{
			template <class Stack>
			struct apply
			{
				typedef stack_iterator<empty_stack> type;
			};
		};

		template <>
		struct clear_impl<stack_tag>
		{
			template <class Seq>
			struct apply : empty_stack {};
		};

		template <>
		struct push_front_impl<stack_tag>
		{
			template <class Stack, class T>
			struct apply : internal::stack_impl<T, Stack> {};
		};

		template <>
		struct pop_front_impl<stack_tag>
		{
			template <class Stack>
			struct apply : Stack::tail {};
		};
	}
}

#endif

