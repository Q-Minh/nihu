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

/** \brief  integrates a kernel over two function spaces and stores the result in a matrix
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

protected:
	/** \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class TestField, class TrialField>
	struct eval_on { struct type {
		typedef TestField test_field_t;		/**< \brief the field type */
		typedef TrialField trial_field_t;		/**< \brief the field type */

		typedef double_integral<kernel_t, test_field_t, trial_field_t> double_integral_t;	/**< \brief the internal integrator type */

		typedef typename double_integral_t::result_t result_t;				/**< \brief result type of field integrator */

		/** \brief evaluate integral on one element type
		 * \param wi constant reference to weighted integral object
		 */
		template <class result_t>
		void operator() (weighted_residual const &wr, result_t &result)
		{
			for (auto test_it = wr.m_test_space.template field_begin<test_field_t>();
				test_it != wr.m_test_space.template field_end<test_field_t>(); ++test_it)
			{
				for (auto trial_it = wr.m_trial_space.template field_begin<trial_field_t>();
					trial_it != wr.m_trial_space.template field_end<trial_field_t>(); ++trial_it)
				{
/*
					result(test_it->get_dofs(),trial_it->get_dofs()) += double_integral_t::eval(*test_it, *trial_it);
*/
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

	/** \brief evaluate integral and return reference to result vector
	 * \return reference to static member result
	 */
	template <class result_t>
	result_t &eval(result_t &result)
	{
		tmp::d_call_each<
			typename test_space_t::field_type_vector_t,
			typename trial_space_t::field_type_vector_t,
			eval_on<tmp::_1, tmp::_2>,
			weighted_residual &,
			result_t &
		>(*this, result);

/*
		typedef field<quad_1_elem, isoparametric_field> test_t;
		typedef field<quad_1_elem, isoparametric_field> trial_t;
		typename eval_on<test_t, trial_t>::type instance;
		instance(*this, result);
*/

		return result;
	}

protected:
	test_space_t const &m_test_space;	/**< \brief reference to the function space */
	trial_space_t const &m_trial_space;	/**< \brief reference to the function space */
};

#endif // ifndef WEIGHTED_RESIDUAL

