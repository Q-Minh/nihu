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
 * \file interval_estimator.hpp
 * \ingroup library
 * \brief speclialisation of complexity_estimator for the interval case
 */
#ifndef INTERVAL_ESTIMATOR_HPP_INCLUDED
#define INTERVAL_ESTIMATOR_HPP_INCLUDED

#include "../core/complexity_estimator.hpp"
#include "../tmp/interval.hpp"

/**
 * \brief specialisation of ::complexity_estimator for the interval case
 * \tparam Interval the intervals to use for the estimation
 */
template <class Interval>
class interval_estimator
{
public:
	/** \brief self returning */
	typedef interval_estimator type;

	template <class test_field_t, class trial_field_t>
	static unsigned eval(
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field
	)
	{
		auto x = test_field.get_elem().get_center();
		auto y = trial_field.get_elem().get_center();
		auto D = trial_field.get_elem().get_linear_size_estimate();
		return eval_interval<Interval>((y - x).norm()/D);
	}
};

template <class Interval1, class Interval2>
struct merge_kernel_complexity_estimators<
	interval_estimator<Interval1>,
	interval_estimator<Interval2>
> : interval_estimator<
	typename merge_intervals<Interval1, Interval2>::type
> {};


#endif // INTERVAL_ESTIMATOR_HPP_INCLUDED
