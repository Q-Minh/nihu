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
 * \file distance_kernel_intervals.hpp
 * \brief Quadrature intervals for distance based complexities
 * \ingroup library
 * \details
 * These intervals are used for selecting the quadrature order for integration
 * based on the non-dimensional distance. The distance can be 
 * non-dimensionalized using the linear size estimate of the elements.
 * 
 * The intervals are specified for different accuracy paramaters, the default 
 * accuracy is given by \ref GLOBAL_ACCURACY. Naturally, higher accuracy takes
 * a higher number of quadrature base points.
 */

#ifndef DISTANCE_KERNEL_INTERVALS_HPP_INCLUDED
#define DISTANCE_KERNEL_INTERVALS_HPP_INCLUDED

#include "../tmp/interval.hpp"
#include "../core/global_definitions.hpp"

namespace NiHu
{

/** 
 * \brief define intervals for distance range and accuracy
 * \tparam Order the reciprocal order
 * \tparam Accuracy -log10(eps) as an integer
 * \details 
 * This is the general base class that can be specialised for different
 * asymptotic orders and accuracy parameters.
 */
template <class asymptotic, unsigned Accuracy = GLOBAL_ACCURACY>
struct distance_kernel_interval;

/** \brief specialisation of ::distance_kernel_interval for log r and 1% error */
template <>
struct distance_kernel_interval<asymptotic::log<1>, 2>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio<13,10>, tmp::integer<int, 4> >,
		tmp::break_point<std::ratio<23,10>, tmp::integer<int, 2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 0> >
	> type;
};


/** \brief specialisation of ::distance_kernel_interval for log r and .1% error */
template <>
struct distance_kernel_interval<asymptotic::log<1>, 3>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio<12,10>, tmp::integer<int, 8> >,
		tmp::break_point<std::ratio<17,10>, tmp::integer<int, 4> >,
		// break_point<std::ratio<51,10>, tmp::integer<int, 2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 2> >
	> type;
};


/** \brief specialisation of ::distance_kernel_interval for 1/r and 1% error */
template <>
struct distance_kernel_interval<asymptotic::inverse<1>, 2>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio<10,10>, tmp::integer<int, 4> >,
		tmp::break_point<std::ratio<28,10>, tmp::integer<int, 2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 0> >
	> type;
};


/** \brief specialisation of ::distance_kernel_interval for 1/r^2 and 1% error */
template <>
struct distance_kernel_interval<asymptotic::inverse<2>, 2>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio<10,10>, tmp::integer<int, 6> >,
		tmp::break_point<std::ratio<15,10>, tmp::integer<int, 4> >,
		tmp::break_point<std::ratio<50,10>, tmp::integer<int, 2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 0> >
	> type;
};


/** \brief specialisation of ::distance_kernel_interval for 1/r^3 and 1% error */
template <>
struct distance_kernel_interval<asymptotic::inverse<3>, 2>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio<10,10>, tmp::integer<int, 8> >,
		tmp::break_point<std::ratio<13,10>, tmp::integer<int, 6> >,
		tmp::break_point<std::ratio<20,10>, tmp::integer<int, 4> >,
		tmp::break_point<std::ratio<70,10>, tmp::integer<int, 2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 0> >
	> type;
};

/** \brief specialisation of ::distance_kernel_interval for 1/r and .1% error */
template <>
struct distance_kernel_interval<asymptotic::inverse<1>, 3>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio<11,10>, tmp::integer<int, 6> >,
		tmp::break_point<std::ratio<18,10>, tmp::integer<int, 4> >,
		tmp::break_point<std::ratio<91,10>, tmp::integer<int, 2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 0> >
	> type;
};

/** \brief specialisation of ::distance_kernel_interval for 1/r^2 and .1% error */
template <>
struct distance_kernel_interval<asymptotic::inverse<2>, 3>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio< 11,10>, tmp::integer<int, 8> >,
		tmp::break_point<std::ratio< 15,10>, tmp::integer<int, 6> >,
		tmp::break_point<std::ratio< 28,10>, tmp::integer<int, 4> >,
		tmp::break_point<std::ratio<158,10>, tmp::integer<int, 2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 0> >
	> type;
};


/** \brief specialisation of ::distance_kernel_interval for 1/r^3 and .1% error */
template <>
struct distance_kernel_interval<asymptotic::inverse<3>, 3>
{
	typedef tmp::vector<
		tmp::break_point<std::ratio< 10,10>, tmp::integer<int,  9> >, // should be 12
		tmp::break_point<std::ratio< 11,10>, tmp::integer<int,  9> >, // should be 10
		tmp::break_point<std::ratio< 14,10>, tmp::integer<int,  8> >,
		tmp::break_point<std::ratio< 19,10>, tmp::integer<int,  6> >,
		tmp::break_point<std::ratio< 36,10>, tmp::integer<int,  4> >,
		tmp::break_point<std::ratio<223,10>, tmp::integer<int,  2> >,
		tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 0> >
	> type;
};

}

#endif /* DISTANCE_KERNEL_INTERVALS_HPP_INCLUDED */
