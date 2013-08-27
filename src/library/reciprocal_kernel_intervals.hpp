#ifndef RECIPROCAL_KERNEL_INTERVALS_HPP_INCLUDED
#define RECIPROCAL_KERNEL_INTERVALS_HPP_INCLUDED

#include "../bem/interval.hpp"

template <unsigned Order, unsigned Error>
struct reciprocal_distance_kernel_interval
{
	// static_assert(false, "reciprocal_distance_kernel_interval is undefined for the selected order and error");
};

template <>
struct reciprocal_distance_kernel_interval<1, 2>
{
	typedef tmp::vector<
		break_point<std::ratio<10,10>, tmp::int_<4> >,
		break_point<std::ratio<28,10>, tmp::int_<2> >,
		break_point<ratio_infinite, tmp::int_<0> >
	> type;
};


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
