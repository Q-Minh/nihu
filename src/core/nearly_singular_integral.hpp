/**
 * \file nearly_singular_integral.hpp 
 * \brief Nearly singular integral general case
 */

#ifndef NIHU_NEARLY_SINGULAR_INTEGRAL_HPP_INCLUDED
#define NIHU_NEARLY_SINGULAR_INTEGRAL_HPP_INCLUDED

#include "field.hpp"
#include "kernel.hpp"

namespace NiHu
{

#if 0
template <class Kernel, class TestField, class TrialField>	
class nearly_singular_distance_needed
{
	static bool needed(
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
};
#endif
	
template <class Kernel, class TestField, class TrialField, class Enable = void>
class nearly_singular_integral
{
public:
	static bool needed(
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
		return false;
	}

	template <class Result>
	static Result &eval(
		Result &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
			return result;
	}
};

} // end of namespace

#endif /* NIHU_NEARLY_SINGULAR_INTEGRAL_HPP_INCLUDED */
