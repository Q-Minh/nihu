#ifndef INTERVAL_ESTIMATOR_HPP_INCLUDED
#define INTERVAL_ESTIMATOR_HPP_INCLUDED

#include "complexity_estimator.hpp"
#include "interval.hpp"

template <class Interval>
class interval_estimator
{
public:
	typedef interval_estimator type;

	template <class test_field_t, class trial_field_t>
	static unsigned eval(
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field
	)
	{
		auto x = test_field.get_elem().get_center();
		auto y = trial_field.get_elem().get_center();
		auto rel_dist = (y-x).norm()/trial_field.get_elem().get_linear_size_estimate();
		return eval_interval<Interval>(rel_dist);
	}
};

template <class Interval1, class Interval2>
struct merge_complexity_estimators<interval_estimator<Interval1>, interval_estimator<Interval2> > :
	interval_estimator<
		typename merge_intervals<Interval1, Interval2>::type
	> {};


#endif // INTERVAL_ESTIMATOR_HPP_INCLUDED
