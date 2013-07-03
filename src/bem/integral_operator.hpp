#ifndef INTEGRAL_OPERATOR_HPP_INCLUDED
#define INTEGRAL_OPERATOR_HPP_INCLUDED

#include "weighted_residual.hpp"
#include "kernel.hpp"

struct local {};
struct non_local {};

template <class Kernel, bool isLocal>
class integral_operator
{
public:
	typedef Kernel kernel_t;
	static bool const is_local = isLocal;

	integral_operator(Kernel const &kernel) :
		m_kernel(kernel)
	{
	}

	Kernel &get_kernel(void)
	{
		return m_kernel;
	}

	template <class TestSpace, class TrialSpace, class Formalism>
	weighted_residual<Formalism, integral_operator, TestSpace, TrialSpace>
		operator()(TestSpace const &test_space, TrialSpace const &trial_space, Formalism)
	{
		return weighted_residual<Formalism, integral_operator, TestSpace, TrialSpace>(*this, test_space, trial_space);
	}

private:
	Kernel m_kernel;
};


template <class Kernel>
integral_operator<Kernel, true> boundary_operator(Kernel const &k, local)
{
	return integral_operator<Kernel, true>(k);
}


template <class Kernel>
integral_operator<Kernel, false> boundary_operator(Kernel const &k, non_local)
{
	return integral_operator<Kernel, false>(k);
}

#endif

