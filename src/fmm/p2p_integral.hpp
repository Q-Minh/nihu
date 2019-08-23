/// \file p2p_integral.hpp

#ifndef FMM_P2P_INTEGRAL_HPP_INCLUDED
#define FMM_P2P_INTEGRAL_HPP_INCLUDED

#include "fmm_operator.hpp"
#include "integral_operator_expression.hpp"
#include "identity_p2p_operator.h"

#include "../core/double_integral.hpp"
#include "../core/single_integral.hpp"

#include "../util/matrix_traits.hpp"
#include "../util/type2tag.hpp"

#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator, class TestField, class TrialField>
class p2p_integral;

template <class Operator, class TestField, class TrialField>
struct integral_operator_expression_traits<p2p_integral<Operator, TestField, TrialField> >
{
	typedef TestField test_input_t;
	typedef TrialField trial_input_t;

	typedef typename NiHu::double_integral<
		typename std::decay<Operator>::type::kernel_t,
		TestField,
		TrialField
	>::result_t result_t;
};


template <class TestField, class TrialField>
struct integral_operator_expression_traits<p2p_integral<identity_p2p_operator, TestField, TrialField> >
{
	typedef TestField test_input_t;
	typedef TrialField trial_input_t;
	typedef typename NiHu::single_integral<
		TestField,
		TrialField
	>::result_t result_t;
};


template <class Operator, class TestField, class TrialField>
class p2p_integral
	: public integral_operator_expression<p2p_integral<Operator, TestField, TrialField> >
	, public fmm_operator<typename std::decay<Operator>::type::fmm_tag>
{
public:
	typedef integral_operator_expression<p2p_integral<Operator, TestField, TrialField> > base_t;

	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	typedef typename std::decay<Operator>::type operator_t;
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	typedef typename operator_t::kernel_t kernel_t;
	typedef NiHu::double_integral<kernel_t, test_field_t, trial_field_t> double_integral_t;

	p2p_integral(Operator &&op, bool sing_check)
		: m_kernel(op.get_kernel())
		, m_singular_check(sing_check)
	{
	}

	size_t rows(test_input_t const &tsi) const
	{
		return num_rows<result_t>::value;
	}

	size_t cols(trial_input_t const &tri) const
	{
		return num_cols<result_t>::value;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		if (m_singular_check)
			return double_integral_t::eval(
				m_kernel, x, y,
				std::true_type());
		else
			return double_integral_t::eval(
				m_kernel, x, y,
				std::false_type());
	}

private:
	kernel_t m_kernel;
	bool const m_singular_check;
};


template <class TestField, class TrialField>
class p2p_integral<identity_p2p_operator, TestField, TrialField>
	: public integral_operator_expression<p2p_integral<identity_p2p_operator, TestField, TrialField> >
	, public fmm_operator<typename identity_p2p_operator::fmm_tag>
{
public:
	typedef integral_operator_expression<p2p_integral<identity_p2p_operator, TestField, TrialField> > base_t;

	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	typedef NiHu::single_integral<test_field_t, trial_field_t> single_integral_t;

	p2p_integral(identity_p2p_operator const &op)
	{
	}

	size_t rows(test_input_t const &tsi) const
	{
		return num_rows<result_t>::value;
	}

	size_t cols(trial_input_t const &tri) const
	{
		return num_cols<result_t>::value;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		/// \todo element coincidence should be checked other way
		bool on_same_elem = &x.get_elem() == &y.get_elem();
		if (on_same_elem)
			return single_integral_t::eval(x, y);
		return result_t::Zero();
	}
};



template <class Operator, class TestTag, class TrialTag>
auto create_p2p_integral(Operator &&op, TestTag, TrialTag, bool sing_check)
{
	typedef typename tag2type<TestTag>::type test_field_t;
	typedef typename tag2type<TrialTag>::type trial_field_t;
	return p2p_integral<Operator, test_field_t, trial_field_t>(std::forward<Operator>(op), sing_check);
}

template <class TestTag, class TrialTag>
auto create_identity_p2p_integral(TestTag, TrialTag)
{
	typedef typename tag2type<TestTag>::type test_field_t;
	typedef typename tag2type<TrialTag>::type trial_field_t;
	return p2p_integral<identity_p2p_operator, test_field_t, trial_field_t>(identity_p2p_operator());
}

} // end of namespace fmm
} // namespace NiHu

#endif
