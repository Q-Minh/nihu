/**
* \file weighted_residual.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class weighted_residual
*/
#ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED
#define WEIGHTED_RESIDUAL_HPP_INCLUDED

#include <Eigen/Dense>

#include <vector>

#include "../tmp/control.hpp"

#include "kernel.hpp"
#include "function_space.hpp"
#include "double_integral.hpp"

/** \brief  integrates a kernel over two function spaces
 * \tparam Kernel the kernel to be integrated
 * \tparam TestSpace the test function space over which integration is performed
 * \tparam TrialSpace the trial function space over which integration is performed
 */
template <class Kernel, class TestSpace, class TrialSpace>
class weighted_residual
{
public:
	typedef TestSpace test_space_t;		/**< \brief template paramter as nested type */
	typedef TrialSpace trial_space_t;	/**< \brief template paramter as nested type */
	typedef Kernel kernel_t;			/**< \brief template paramter as nested type */

private:
	/** \brief evaluate weighted residual on homogeneous function spaces
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam TestField the test field type over which integration is performed
	 * \tparam TrialField the trial field type over which integration is performed
	 */
	template <class TestField, class TrialField>
	struct eval_on { struct type {
		typedef TestField test_field_t;		/**< \brief the test field type */
		typedef TrialField trial_field_t;	/**< \brief the trial field type */

		typedef double_integral<kernel_t, test_field_t, trial_field_t> double_integral_t;	/**< \brief the internal integrator type */
		typedef typename double_integral_t::result_t result_t;				/**< \brief result type of field integrator */

		template <class result_t>
		void operator() (weighted_residual const &wr, result_t &result)
		{
			for (auto test_it = wr.m_test_space.template field_begin<test_field_t>();
				test_it != wr.m_test_space.template field_end<test_field_t>(); ++test_it)
			{
				for (auto trial_it = wr.m_trial_space.template field_begin<trial_field_t>();
					trial_it != wr.m_trial_space.template field_end<trial_field_t>(); ++trial_it)
				{
					// result(test_it->get_dofs(),trial_it->get_dofs()) += double_integral_t::eval(*test_it, *trial_it);
					std::cout << double_integral_t::eval(*test_it, *trial_it) << std::endl << std::endl;
				}
			}
		}
	};};

public:
	/** \brief constructor initialises the function space reference members
	 * \param test_space the function space over which integration is performed
	 * \param trial_space the function space over which integration is performed
	 */
	weighted_residual(test_space_t const &test_space, trial_space_t const &trial_space) : m_test_space(test_space), m_trial_space(trial_space)
	{
	}

	/** \brief evaluate weighted residual and return reference to the result matrix
	 * \details This function evaluates the weighted residual of a kernel over a test and a trial field.
	 * The two fields have been initialised by the constructor, and are stored in member variables.
	 * The result is copied into a user specified structure, passed as function argument.
	 *
	 * \tparam result_t type of the result matrix
	 * \param [out] result reference to the result matrix
	 * \return reference to the result matrix for cascading
	 */
	template <class result_t>
	result_t &eval(result_t &result)
	{
		// Integration is performed separately on homogneous subfields, using tmp::d_call_each
		// d_call_each calls eval_on for each element of the the Descartes product of the test and
		// trial field type vectors
		tmp::d_call_each<
			typename test_space_t::field_type_vector_t,
			typename trial_space_t::field_type_vector_t,
			eval_on<tmp::_1, tmp::_2>,
			weighted_residual const &,
			result_t &
		>(*this, result);

		return result;
	}

protected:
	test_space_t const &m_test_space;	/**< \brief reference to the function space */
	trial_space_t const &m_trial_space;	/**< \brief reference to the function space */
};

#endif // ifndef WEIGHTED_RESIDUAL

