/** \file integral_operator.hpp
* \brief declaration of class ::integral_operator
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/

#ifndef INTEGRAL_OPERATOR_HPP_INCLUDED
#define INTEGRAL_OPERATOR_HPP_INCLUDED

#include "dirac_wrapper.hpp"
#include "weighted_residual.hpp"

/** \brief collection of options representing integral operator behavior */
namespace operator_option
{
	/** \brief the operator acts within one element, resulting in sparse matrices */
	struct local {};
	/** \brief the operator acts between different elements, resulting in dense matrices */
	struct non_local {};
}

/** \brief Integral operator representation
* \tparam Kernel the kernel class the operator evaluates
* \tparam LocalOption indicates if the operator operates locally or non_locally
*/
template <class Kernel, class LocalOption>
class integral_operator
{
public:
	/** \brief template argument as nested type */
	typedef Kernel kernel_t;
	/** \brief indicates if the operator operates locally or not */
	static bool const is_local = std::is_same<LocalOption, operator_option::local>::value;

	/** \brief constructor from kernel reference
	* \param [in] kernel reference to the kernel
	*/
	integral_operator(Kernel const &kernel) :
		m_kernel(kernel)
	{
	}

	/** \brief return kernel reference
	* \return reference to the kernel
	*/
	Kernel &get_kernel(void)
	{
		return m_kernel;
	}

	/** \brief evaluate operator on a test and a trial function space
	* \tparam TestSpace type of the test function space
	* \tparam TrialSpace type of the trial function space
	* \param test_space the test function space
	* \param trial_space the trial function space
	*/
	template <class TestSpace, class TrialSpace>
	weighted_residual<formalism::general, integral_operator, TestSpace, TrialSpace>
		operator()(TestSpace const &test_space, TrialSpace const &trial_space)
	{
		return weighted_residual<formalism::general, integral_operator, TestSpace, TrialSpace>(
			*this, test_space, trial_space);
	}


	/** \brief evaluate operator on a test and a trial function space
	* \tparam TestSpace type of the test function space
	* \tparam TrialSpace type of the trial function space
	* \param test_wrapper the Dirac wrapper of the test function space
	* \param trial_space the trial function space
	*/
	template <class TestSpace, class TrialSpace>
	weighted_residual<formalism::collocational, integral_operator, TestSpace, TrialSpace>
		operator()(dirac_wrapper<TestSpace> const &test_wrapper, TrialSpace const &trial_space)
	{
		return weighted_residual<formalism::collocational, integral_operator, TestSpace, TrialSpace>(
			*this, test_wrapper.get_function_space(), trial_space);
	}


	/** \brief evaluate operator on a pair of Dirac function spaces
	* \tparam TestSpace type of the test function space
	* \param test_wrapper the dirac wrapper of the function space
	*/
	template <class TestSpace>
	weighted_residual<formalism::full_dirac, integral_operator, TestSpace, TestSpace>
		operator()(dirac_wrapper<TestSpace> const &test_wrapper)
	{
		return weighted_residual<formalism::full_dirac, integral_operator, TestSpace, TestSpace>(
			*this, test_wrapper.get_function_space(), test_wrapper.get_function_space());
	}


private:
	/** \brief he underlying kernel */
	Kernel m_kernel;
};


/** \brief factory function of an integral operator
* \tparam Kernel the kernel type
* \tparam LocalOption the kernel operation option
* \param [in] kernel the kernel
* \return the integral operator object
*/
template <class Kernel, class LocalOption>
integral_operator<Kernel, LocalOption>
	create_integral_operator(Kernel const &kernel, LocalOption)
{
	return integral_operator<Kernel, LocalOption>(kernel);
}


/** \brief factory function of a (default) non_local integral operator
* \tparam Kernel the kernel type
* \param [in] kernel the kernel
* \return the integral operator object
*/
template <class Kernel>
integral_operator<Kernel, operator_option::non_local>
	create_integral_operator(Kernel const &kernel)
{
	return integral_operator<Kernel, operator_option::non_local>(kernel);
}

#endif

