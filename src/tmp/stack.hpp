// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

	/** \brief the stack iterator type used by tmp::begin, tmp::end and tmp::deref */
	template <class Seq>
	struct stack_iterator
	{
		/** \brief self-returning matafunction */
		typedef stack_iterator type;
	};

	/** \brief increment a stack iterator */
	template <class Seq>
	struct next<stack_iterator<Seq> > : stack_iterator<typename Seq::tail> {};

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
			struct apply : stack_iterator<Stack> {};
		};

		template <>
		struct end_impl<stack_tag>
		{
			template <class Stack>
			struct apply : stack_iterator<empty_stack> {};
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

