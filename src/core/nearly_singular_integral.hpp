#ifndef NIHU_NEARLY_SINGULAR_INTEGRAL_SHORTCUT_HPP_INCLUDED
#define NIHU_NEARLY_SINGULAR_INTEGRAL_SHORTCUT_HPP_INCLUDED

#include "field.hpp"
#include "kernel.hpp"

namespace NiHu
{

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

}

#endif // NIHU_NEARLY_SINGULAR_INTEGRAL_SHORTCUT_HPP_INCLUDED
