/**
* \file weighted_residual.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class ::weighted_residual
*/
#ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED
#define WEIGHTED_RESIDUAL_HPP_INCLUDED

#include "matrix_block.hpp"
#include "function_space.hpp"
#include "double_integral.hpp"
#include "../util/dual_range.hpp"

/** \brief  integrates a kernel over two function spaces
 * \tparam Kernel the kernel to be integrated
 * \tparam TestSpace the test function space over which integration is performed
 * \tparam TrialSpace the trial function space over which integration is performed
 */
template <class Formalism, class Operator, class TestSpace, class TrialSpace>
class weighted_residual
{
public:
	/** \brief template parameter as nested type */
	typedef Operator operator_t;
	/** \brief template parameter as nested type */
	typedef TestSpace test_space_t;
	/** \brief template parameter as nested type */
	typedef TrialSpace trial_space_t;

	/** \brief template parameter as nested type */
	typedef typename operator_t::kernel_t kernel_t;
	/** \brief indicates if operator is local or not */
	static bool const is_local = operator_t::is_local;

private:

	/** \brief evaluate weighted residual on homogeneous function spaces
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam TestField the test field type over which integration is performed
	 * \tparam TrialField the trial field type over which integration is performed
	 */
	template <class TestField, class TrialField>
	struct eval_on { struct type {
		/** \brief the test field type */
		typedef TestField test_field_t;
		/** \brief the trial field type */
		typedef TrialField trial_field_t;

		/** \brief indicates if the two field types are defined over the same element type */
		static bool const same_elem_type = std::is_same<
			typename test_field_t::elem_t,
			typename trial_field_t::elem_t
		>::value;

		/** \brief the internal integrator type */
		typedef double_integral<kernel_t, test_field_t, trial_field_t> double_integral_t;
		/** \brief result type of field integrator */
		typedef typename double_integral_t::result_t result_t;

		/** \brief evaluate weighted residual on homogeneous function spaces
		 * \tparam result_t the type of the result matrix
		 * \param [out] result reference to the result matrix the result block is inserted to
		 * \param [in] kernel the kernel to integrate
		 * \param [in] test_space the test function space
		 * \param [in] trial_space the trial function space
		 */
		template <class result_t>
		void operator() (result_t &result,
			kernel_t &kernel,
			test_space_t const &test_space,
			trial_space_t const &trial_space)
		{
			auto rng = create_dual_range(
				iteration::plain(),
				test_space.template field_begin<test_field_t>(),
				test_space.template field_end<test_field_t>(),
				trial_space.template field_begin<trial_field_t>(),
				trial_space.template field_end<trial_field_t>()
			);

			auto it = rng.begin();
			auto end = rng.end();
			while (it != end)
			{
				if (!is_local ||
					( same_elem_type && (*it.outer()).get_elem().get_id() == (*it.inner()).get_elem().get_id()) )
				{
					block(result, (*it.outer()).get_dofs(), (*it.inner()).get_dofs())
						+= double_integral_t::eval(Formalism(),
							kernel, *it.outer(), *it.inner());
				}
				++it;
			}
		}
	};};

public:
	/** \brief constructor from operator and spaces
	* \param [in] op the integral operator
	* \param test_space the test function space
	* \param trial_space the trial function space
	*/
	weighted_residual(
		operator_t &op,
		test_space_t const &test_space,
		trial_space_t const &trial_space) :
		m_operator(op), m_test_space(test_space), m_trial_space(trial_space)
	{
	}

	/** \brief evaluate weighted residual and return reference to the result matrix
	 * \tparam result_t type of the result matrix
	 * \param [out] result reference to the result matrix
	 * \return reference to the result matrix for cascading
	 */
	template <class result_t>
	result_t &eval(result_t &result)
	{
		tmp::d_call_each<
			typename test_space_t::field_type_vector_t,
			typename trial_space_t::field_type_vector_t,
			eval_on<tmp::_1, tmp::_2>,
			result_t &,
			kernel_t &,
			test_space_t const &,
			trial_space_t const &
		>(result, m_operator.get_kernel(), m_test_space, m_trial_space);

		return result;
	}

private:
	/** \brief the integral operator reference */
	operator_t &m_operator;
	/** \brief the test space reference */
	test_space_t const &m_test_space;
	/** \brief the trial space reference */
	trial_space_t const &m_trial_space;
};


template <class Result, class WR>
void eval_into(Result &res, WR wr)
{
	wr.eval(res);
}


#endif // ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED

