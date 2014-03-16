// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
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
* \file control.hpp
* \ingroup tmp
* \brief Implementation of code generating control structures
*/
#ifndef CONTROL_HPP_INCLUDED
#define CONTROL_HPP_INCLUDED

#include <type_traits>

#include "sequence.hpp"
#include "lambda.hpp"
#include "algorithm.hpp"

/** \brief template metaprogramming functions */
namespace tmp
{
	/** \brief implementation namespace */
	namespace internal
	{
		/**
		* \brief implementation of ::call_each
		* \tparam Begin the begin iterator of the sequence
		* \tparam End the end iterator of the sequence
		* \tparam Transform a metafunction with a callable internal struct called type
		* \tparam Args function arguments
		*/
		template <class Begin, class End, class Transform, class...Args>
		struct call_each_impl
		{
			static void eval(Args ...args)
			{
				// instantiate and call actual templated version of Transform
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				c(args...);
				// envoke call_each for the remaining part of the sequence
				call_each_impl<
					typename next<Begin>::type,
					End,
					Transform,
					Args...
				>::eval(args...);
			}
		};

		/** \brief specialisation of ::call_each_impl for the terminating case */
		template <class End, class Transform, class...Args>
		struct call_each_impl<End, End, Transform, Args...>
		{
			static void eval(Args...) { }
		};
	}

	/**
	* \brief Call Transform<x>::type(Args...) for each element x of the sequence Seq
	* \tparam Seq the sequence of values x
	* \tparam Transform a metafunction with a callable internal struct called type
	* \tparam Args function arguments
	*/
	template <class Seq, class Transform, class...Args>
	static void call_each(Args...args)
	{
		internal::call_each_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform,
			Args...
		>::eval(args...);
	}

	namespace internal
	{
		/**
		* \brief implementation of d_call_each
		* \tparam Begin begin iterator of the outer sequence
		* \tparam End end iterator of the outer sequence
		* \tparam SeqIn inner type sequence
		* \tparam Transform binary metafunction with a callable internal struct called type
		* \tparam Args function arguments
		*/
		template <class Begin, class End, class SeqIn, class Transform, class...Args>
		struct d_call_each_impl
		{
			// bind transform to the first argument
			typedef typename lambda<Transform>::type::template apply<
				typename deref<Begin>::type, _1
			> partially_evaluated_transform;

			static void eval(Args...args)
			{
				// call call_each with partially evaluated Transform
				call_each<
					SeqIn,
					partially_evaluated_transform,
					Args...
				>(args...);
				// call d_call_each with the rest of the outer sequence
				d_call_each_impl<
					typename next<Begin>::type,
					End,
					SeqIn,
					Transform,
					Args...
				>::eval(args...);
			}
		};


		/** \brief terminating case of d_call_each_impl */
		template <class End, class SeqIn, class Transform, class...Args>
		struct d_call_each_impl<End, End, SeqIn, Transform, Args...>
		{
			static void eval(Args...) {}
		};
	}


	/**
	* \brief Call Transform<x,y>::type(Args...) for each element x,y of the Descartes product of two sequences
	* \tparam SeqOut outer type sequence
	* \tparam SeqIn inner type sequence
	* \tparam Transform binary metafunction with a callable internal struct called type
	* \tparam Args function arguments
	*/
	template <class SeqOut, class SeqIn, class Transform, class...Args>
	static void d_call_each(Args...args)
	{
		internal::d_call_each_impl<
			typename begin<SeqOut>::type,
			typename end<SeqOut>::type,
			SeqIn,
			Transform,
			Args...
		>::eval(args...);
	}

	namespace internal
	{
		/**
		* \brief implementation of ::call_until
		* \tparam Begin the begin iterator of the sequence
		* \tparam End the end iterator of the sequence
		* \tparam Transform a metafunction with a callable internal struct called type
		* \tparam Args function arguments
		*/
		template <class Begin, class End, class Transform, class...Args>
		struct call_until_impl
		{
			static bool eval(Args...args)
			{
				typedef typename apply<Transform, typename deref<Begin>::type>::type cur;
				cur c;
				if (!c(args...))
				{
					return call_until_impl<
						typename next<Begin>::type,
						End,
						Transform,
						Args...
					>::eval(args...);
				}
				return true;
			}
		};

		/** \brief terminating case of call_until_impl */
		template <class End, class Transform, class...Args>
		struct call_until_impl<End, End, Transform, Args...>
		{
			static bool eval(Args...) { return false; }
		};
	}

	/**
	 * \brief Call Transform<x>::type(Args...) for each element x of the sequence Seq until it returns true
	 * \tparam Seq the sequence of values x
	 * \tparam Transform a metafunction with a callable internal struct called type
	 * \tparam Args function arguments
	 */
	template <class Seq, class Transform, class...Args>
	static bool call_until(Args...args)
	{
		return internal::call_until_impl<
			typename begin<Seq>::type,
			typename end<Seq>::type,
			Transform,
			Args...
		>::eval(args...);
	}
}

#endif // CONTROL_HPP_INCLUDED

