#ifndef COMPLEXITY_ESTIMATOR_HPP_INCLUDED
#define COMPLEXITY_ESTIMATOR_HPP_INCLUDED

#include "field.hpp"
#include "../tmp/integer.hpp"

class interval_estimator
{
public:
	template <class test_field_t, class trial_field_t>
	static unsigned eval(
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field
	)
	{
		return 2;
	}
};

template <class TestField, class TrialField, class KernelEstimator>
class complexity_estimator
{
public:
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	enum {
		test_field_complexity =
			shape_set_traits<typename test_field_t::nset_t>::polynomial_order +
			shape_set_traits<typename test_field_t::lset_t>::jacobian_order,
		trial_field_complexity =
			shape_set_traits<typename trial_field_t::nset_t>::polynomial_order +
			shape_set_traits<typename trial_field_t::lset_t>::jacobian_order
	};

	static unsigned const total_field_complexity = tmp::max_<
		std::integral_constant<unsigned, test_field_complexity>,
		std::integral_constant<unsigned, trial_field_complexity>
	>::value;

	static unsigned eval(
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field
	)
	{
		return total_field_complexity +
			KernelEstimator::eval(test_field, trial_field);
	}
};


template <class Estim1, class...Estimators>
struct merge_kernel_complexity_estimators
{
	struct type
	{
		template <class test_field_t, class trial_field_t>
		static unsigned eval(
			field_base<test_field_t> const &test_field,
			field_base<trial_field_t> const &trial_field
		)
		{
			return std::max(
				Estim1::eval(test_field, trial_field),
				merge_kernel_complexity_estimators<Estimators...>::type::eval(
					test_field,
					trial_field
				)
			);
		}	// end of function
	};	// end of struct type
};

template <class Estim1>
struct merge_kernel_complexity_estimators<Estim1>
{
	struct type
	{
		template <class test_field_t, class trial_field_t>
		static unsigned eval(
			field_base<test_field_t> const &test_field,
			field_base<trial_field_t> const &trial_field
		)
		{
			return Estim1::eval(test_field, trial_field);
		}
	};
};

#endif
