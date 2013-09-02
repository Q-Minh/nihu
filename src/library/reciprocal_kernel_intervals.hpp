/**
 * \file reciprocal_kernel_intervals.hpp
 * \brief quadrature intervals for the 1/r^(order) type kernels
 */

#ifndef RECIPROCAL_KERNEL_INTERVALS_HPP_INCLUDED
#define RECIPROCAL_KERNEL_INTERVALS_HPP_INCLUDED

#include "../core/interval.hpp"

/** \brief define intervals for reciprocal order and accuracy
 * \tparam Order the reciprocal order
 * \tparam Accuracy -log10(eps) as an integer
 */
template <unsigned Order, unsigned Accuracy>
struct reciprocal_distance_kernel_interval;


/** \brief specialisation of ::reciprocal_distance_kernel_interval for 1/r and 1% error */
template <>
struct reciprocal_distance_kernel_interval<1, 2>
{
	typedef tmp::vector<
		break_point<std::ratio<10,10>, tmp::int_<4> >,
		break_point<std::ratio<28,10>, tmp::int_<2> >,
		break_point<ratio_infinite, tmp::int_<0> >
	> type;
};


/** \brief specialisation of ::reciprocal_distance_kernel_interval for 1/r^2 and 1% error */
template <>
struct reciprocal_distance_kernel_interval<2, 2>
{
	typedef tmp::vector<
		break_point<std::ratio<10,10>, tmp::int_<6> >,
		break_point<std::ratio<15,10>, tmp::int_<4> >,
		break_point<std::ratio<50,10>, tmp::int_<2> >,
		break_point<ratio_infinite, tmp::int_<0> >
	> type;
};


/** \brief specialisation of ::reciprocal_distance_kernel_interval for 1/r^3 and 1% error */
template <>
struct reciprocal_distance_kernel_interval<3, 2>
{
	typedef tmp::vector<
		break_point<std::ratio<10,10>, tmp::int_<8> >,
		break_point<std::ratio<13,10>, tmp::int_<6> >,
		break_point<std::ratio<20,10>, tmp::int_<4> >,
		break_point<std::ratio<70,10>, tmp::int_<2> >,
		break_point<ratio_infinite, tmp::int_<0> >
	> type;
};




/** \brief specialisation of ::reciprocal_distance_kernel_interval for 1/r and .1% error */
template <>
struct reciprocal_distance_kernel_interval<1, 3>
{
	typedef tmp::vector<
		break_point<std::ratio<11,10>, tmp::int_<6> >,
		break_point<std::ratio<18,10>, tmp::int_<4> >,
		break_point<std::ratio<91,10>, tmp::int_<2> >,
		break_point<ratio_infinite, tmp::int_<0> >
	> type;
};


/** \brief specialisation of ::reciprocal_distance_kernel_interval for 1/r^2 and .1% error */
template <>
struct reciprocal_distance_kernel_interval<2, 3>
{
	typedef tmp::vector<
		break_point<std::ratio< 11,10>, tmp::int_<8> >,
		break_point<std::ratio< 15,10>, tmp::int_<6> >,
		break_point<std::ratio< 28,10>, tmp::int_<4> >,
		break_point<std::ratio<158,10>, tmp::int_<2> >,
		break_point<ratio_infinite, tmp::int_<0> >
	> type;
};


/** \brief specialisation of ::reciprocal_distance_kernel_interval for 1/r^3 and .1% error */
template <>
struct reciprocal_distance_kernel_interval<3, 3>
{
	typedef tmp::vector<
		break_point<std::ratio< 10,10>, tmp::int_< 9> >, // should be 12
		break_point<std::ratio< 11,10>, tmp::int_< 9> >, // should be 10
		break_point<std::ratio< 14,10>, tmp::int_< 8> >,
		break_point<std::ratio< 19,10>, tmp::int_< 6> >,
		break_point<std::ratio< 36,10>, tmp::int_< 4> >,
		break_point<std::ratio<223,10>, tmp::int_< 2> >,
		break_point<ratio_infinite, tmp::int_<0> >
	> type;
};


#endif // RECIPROCAL_KERNEL_INTERVALS_HPP_INCLUDED
