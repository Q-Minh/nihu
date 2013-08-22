#ifndef INTERVAL_HPP_INCLUDED
#define INTERVAL_HPP_INCLUDED

#include "tmp/vector.hpp"
#include <ratio>
#include "tmp/integer.hpp"

template <class X, class Y>
struct break_point
{
	/** \brief self-returning */
	typedef break_point type;
	/** \brief the x data */
	typedef X x;
	/** \brief the y data */
	typedef Y y;

	/** \brief convert rational x value to scalar(double) */
	template<class scalar_t = double>
	constexpr static scalar_t x_value(void)
	{
		return scalar_t(x::num)/scalar_t(x::den);
	}
};

template <class BreakPoint>
struct get_x;

template <class X, class Y>
struct get_x<break_point<X, Y> >
{
	typedef X type;
};

template <class InterPoint>
struct get_y;

template <class X, class Y>
struct get_y<break_point<X, Y> >
{
	typedef Y type;
};




#endif // INTERVAL_HPP_INCLUDED
