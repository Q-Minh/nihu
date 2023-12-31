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
 * \file interval.hpp
 * \brief A compile time interval
 * \ingroup tmp
 * \details
 * Intervals are used for evaluating qudrature orders when integrating 
 * distance dependent kernels.
 */
#ifndef INTERVAL_HPP_INCLUDED
#define INTERVAL_HPP_INCLUDED

#include "../core/global_definitions.hpp"

#include "algorithm.hpp"
#include "integer.hpp"
#include "vector.hpp"

#include <ratio>
#include <stdexcept>

namespace tmp
{

/** \brief define giga as the infinite as it fits into 32 bits */
typedef std::giga ratio_infinite;

/**
 * \brief A break point consisting of a X and a Y value
 * \tparam X the x value type
 * \tparam Y the y value type
 * \details 
 * A break point can represent a non-dimensionalized distance limit 
 * where the number of quadrature base points needed for the  precise evalution 
 * of the integral changes. A list of break points make up the interval. The x 
 * value of the break point is non-dimensional distance limit, and the y value 
 * is the number of base points required.
 */
template <class X, class Y>
struct break_point
{
	/** \brief Self-returning */
	typedef break_point type;
	/** \brief the x data */
	typedef X x;
	/** \brief the y data */
	typedef Y y;

	/**
	 * \brief Convert rational x value to scalar(double)
	 * \tparam scalar_t The conversion type (double is default)
	 * \return The rational number converted to scalar_t
	 */
	template<class scalar_t = double>
	constexpr static scalar_t x_value(void)
	{
		return scalar_t(x::num)/scalar_t(x::den);
	}
};


/**
 * \brief Merge two intervals
 * \tparam Inter1 the first interval
 * \tparam Inter2 the second interval
 * \details 
 * The result of the merge is a new interval containing all the breakpoints
 * of both intervals. The x values are unique-sorted. If an x value appears
 * in both intervals, the merged interval will contain the higher y value
 * for that x value.
 */
template <class Inter1, class Inter2>
struct merge_intervals
{
	/**
	 * \brief Copy condition when zipping an interval
	 * \tparam Iter The iterator pointing to the source element
	 * \tparam Begin Tterator pointing to the first element
	 */
	template <class Iter, class Begin>
	struct copy_cond
	{
		/** \brief previous element */
		typedef typename tmp::deref<typename tmp::prev<Iter>::type>::type prev;
		/** \brief current element */
		typedef typename tmp::deref<Iter>::type cur;
		/** \brief copy if cur.x > prev.x and cur.y < prev.y */
		typedef typename std::integral_constant<
			bool,
			std::ratio_greater<
				typename cur::x,
				typename prev::x
			>::type::value
			&&
			tmp::less<
				typename cur::y,
				typename prev::y
			>::type::value
		> type;
	};

	/** 
	 * \brief Copy condition when zipping an interval
	 * \details 
	 * Specialisation for the first element (always copy)
	 */
	template <class Begin>
	struct copy_cond<Begin, Begin> : std::true_type {};

	/** \brief comparison condition to sort in descending order by y */
	template <class BP1, class BP2>
	struct compare_by_y_desc : tmp::less<
		typename BP1::y,
		typename BP2::y
	> {};

	/** \brief comparison condition to sort in ascending order by x */
	template <class BP1, class BP2>
	struct compare_by_x_asc : std::ratio_less<
		typename BP1::x,
		typename BP2::x
	> {};

	/** \brief sort the concatenated intervals by x ascending */
	typedef typename tmp::bubble_sort<
		typename tmp::concatenate<Inter1, Inter2>::type,
		compare_by_x_asc<tmp::_1, tmp::_2>
	>::type half_sorted;

	/** \brief sort by y descending */
	typedef typename tmp::bubble_sort<
		half_sorted,
		compare_by_y_desc<tmp::_1, tmp::_2>
	>::type sorted;

	/** \brief zip using the copy condition */
	typedef typename tmp::copy_if<
		sorted,
		tmp::inserter<
			typename tmp::empty<sorted>::type,
			tmp::push_back<tmp::_1, tmp::_2>
		>,
        copy_cond<tmp::_1, typename tmp::begin<sorted>::type>
	>::type type;
};

/** 
 * \brief Evaluate an interval at a given distance
 * \tparam interval the interval to evaluate
 * \param [in] r the scalar distance
 * \return the appropriate interval limit
 * \details 
 * Returns the y value of the first breakpoint of the interval whose x value
 * is greater than the non-dimensional distance \c r.
 */
template <class interval>
int eval_interval(double r)
{
	/** \brief get first break point */
	typedef typename tmp::deref<
		typename tmp::begin<interval>::type
	>::type bp;
	// compare radius with first x value
	if (r < bp::x_value())
		return bp::y::value;
	else
		return eval_interval<typename tmp::pop_front<interval>::type>(r);
}

/** \brief error terminating case of eval_interval */
template <>
int eval_interval<tmp::vector<> >(double);

} // end of namespace tmp

#endif /* INTERVAL_HPP_INCLUDED */

