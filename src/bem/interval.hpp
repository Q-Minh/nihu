#ifndef INTERVAL_HPP_INCLUDED
#define INTERVAL_HPP_INCLUDED

#include "tmp/integer.hpp"
#include <ratio>
#include "tmp/vector.hpp"
#include "tmp/algorithm.hpp"

/** \brief define exa as the infinite */
typedef std::exa ratio_infinite;

/**
 * \brief a break point consisting of a X and a Y value
 * \tparam X the x value type
 * \tparam Y the y value type
 */
template <class X, class Y>
struct break_point
{
	/** \brief self-returning */
	typedef break_point type;
	/** \brief the x data */
	typedef X x;
	/** \brief the y data */
	typedef Y y;

	/**
	 * \brief convert rational x value to scalar(double)
	 * \tparam scalar_t the conversion type (double is default)
	 * \return the rational number converted to scalar_t
	 */
	template<class scalar_t = double>
	constexpr static scalar_t x_value(void)
	{
		return scalar_t(x::num)/scalar_t(x::den);
	}
};


/**
 * \brief merge two intervals
 * \tparam Inter1 the first interval
 * \tparam Inter2 the second interval
 */
template <class Inter1, class Inter2>
struct merge_intervals
{
	/**
	 * \brief copy condition when zipping an interval
	 * \tparam Iter the iterator pointing to the source element
	 * \tparam Begin iterator pointing to the first element
	 */
	template <class Iter, class Begin>
	struct copy_cond
	{
		// previous and current elements
		typedef typename tmp::deref<typename tmp::prev<Iter>::type>::type prev;
		typedef typename tmp::deref<Iter>::type cur;
		// copy if cur.x > prev.x and cur.y < prev.y
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

	/** \brief specialisation for the first element */
	template <class Begin>
	struct copy_cond<Begin, Begin> : std::true_type {};

	/** \brief comparison condition to sort in descending order by y */
	template <class BP1, class BP2>
	struct compare_by_y_desc : tmp::less<
		typename BP1::y,
		typename BP2::y
	> {};

	/** \brief comparison condition to sort in descending order by x */
	template <class BP1, class BP2>
	struct compare_by_x_asc : std::ratio_less<
		typename BP1::x,
		typename BP2::x
	> {};

	// sort the concatenated intervals by x ascending
	typedef typename tmp::bubble_sort<
		typename tmp::concatenate<Inter1, Inter2>::type,
		compare_by_x_asc<tmp::_1, tmp::_2>
	>::type half_sorted;

	// sort by y descending
	typedef typename tmp::bubble_sort<
		half_sorted,
		compare_by_y_desc<tmp::_1, tmp::_2>
	>::type sorted;

	// zip using the copy condition
	typedef typename tmp::copy_if<
		sorted,
		tmp::inserter<
			typename tmp::empty<sorted>::type,
			tmp::push_back<tmp::_1, tmp::_2>
		>,
        copy_cond<tmp::_1, typename tmp::begin<sorted>::type>
	>::type type;
};

template <class interval>
int eval_interval(double r)
{
	// get first break point
	typedef typename tmp::deref<
		typename tmp::begin<interval>::type
	>::type bp;
	// compare radius with first x value
	if (r < bp::x_value())
		return bp::y::value;
	else
		return eval_interval<typename tmp::pop_front<interval>::type>(r);
}

template <>
int eval_interval<tmp::vector<> >(double)
{
	throw "cannot determine interval value";
}


#endif // INTERVAL_HPP_INCLUDED
