#ifndef NORMAL_DERIVATIVE_SINGULAR_INTEGRALS_HPP_INCLUDED
#define NORMAL_DERIVATIVE_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "normal_derivative_kernel.hpp"
#include "../core/singular_integral_shortcut.hpp"

namespace NiHu
{

template <class DK, class TestField, class TrialField>
class singular_integral_shortcut<
	normal_derivative_kernel<DK, 0, 1>,
	TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static constexpr result_t &eval(
		result_t &result,
		kernel_base<normal_derivative_kernel<DK, 0, 1> > const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};


template <class DK, class TestField, class TrialField>
class singular_integral_shortcut<
	normal_derivative_kernel<DK, 1, 0>,
	TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static constexpr result_t &eval(
		result_t &result,
		kernel_base<normal_derivative_kernel<DK, 1, 0> > const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};


template <class DK, class TestField, class TrialField>
class singular_integral_shortcut<
	normal_derivative_kernel<DK, 0, 1>,
	TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, line_1_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static constexpr result_t &eval(
		result_t &result,
		kernel_base<normal_derivative_kernel<DK, 0, 1> > const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};


template <class DK, class TestField, class TrialField>
class singular_integral_shortcut<
	normal_derivative_kernel<DK, 1, 0>,
	TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, line_1_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static constexpr result_t &eval(
		result_t &result,
		kernel_base<normal_derivative_kernel<DK, 0, 1> > const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};


} // end of namespace NiHu

#endif // NORMAL_DERIVATIVE_SINGULAR_INTEGRALS_HPP_INCLUDED
