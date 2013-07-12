/**
 * \file assembly.hpp
 * \brief assemble WR matrices from field wr submatrices
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef ASSEMBLY_HPP_INCLUDED
#define ASSEMBLY_HPP_INCLUDED

#include "../util/matrix_block.hpp"
#include "integral_operator.hpp"

template <class TestSpace, class TrialSpace>
class assembly
{
private:
	/** \brief implementation of ::eval_on
	 * \tparam isEmpty indicates if the operation is void or not
	 * \tparam Operator the integral operator type that is applied
	 * \tparam TestField the test field type over which integration is performed
	 * \tparam TrialField the trial field type over which integration is performed
	 */
	template <bool is_empty, class Operator, class TestField, class TrialField>
	struct eval_on_impl
	{ struct type {
		/** \brief evaluate weighted residual on homogeneous function spaces
		 * \tparam result_t the type of the result matrix
		 * \param [out] result reference to the result matrix the result block is inserted to
		 * \param [in] op the integral operator to apply
		 * \param [in] test_space the test function space
		 * \param [in] trial_space the trial function space
		 */
		template <class result_t>
		void operator() (result_t &result,
			Operator const &op,
			function_space_base<TestSpace> const &test_space,
			function_space_base<TrialSpace> const &trial_space)
		{
			for (auto test_it = test_space.template field_begin<TestField>();
				test_it != test_space.template field_end<TestField>(); ++test_it)
				for (auto trial_it = trial_space.template field_begin<TrialField>();
					trial_it != trial_space.template field_end<TrialField>(); ++trial_it)
					block(result, test_it->get_dofs(), trial_it->get_dofs())
						+= op.eval_on_fields(*test_it, *trial_it);
		}
	};};


	/** \brief trivial specialisation of ::eval_on_impl for the empty case */
	template <class Operator, class TestField, class TrialField>
	struct eval_on_impl<true, Operator, TestField, TrialField>
	{ struct type {
		template <class result_t>
		void operator() (
			result_t &,
			Operator const &,
			function_space_base<TestSpace> const &,
			function_space_base<TrialSpace> const &)
		{
		}
	};};

	/** \brief evaluate weighted residual on homogeneous function spaces
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam Operator the integral operator type that is applied
	 * \tparam TestField the test field type over which integration is performed
	 * \tparam TrialField the trial field type over which integration is performed
	 */
	template <class Operator, class TestField, class TrialField>
	struct eval_on : eval_on_impl<
		integral_operator_traits<Operator>::is_local && !std::is_same<typename TestField::elem_t, typename TrialField::elem_t>::value,
		Operator, TestField, TrialField
	> {};


public:
	/** \brief evaluate weighted residual into result matrix
	* \tparam result_t the result matrix type
	* \tparam Operator the operator type
	* \param [out] result the result matrix reference
	* \param [in] op reference to the integral operator
	* \param [in] test reference to the test function space
	* \param [in] trial reference to the trial function space
	*/
	template <class result_t, class Operator>
	static result_t &eval_into(
		result_t &result,
		Operator const &op,
		function_space_base<TestSpace> const &test,
		function_space_base<TrialSpace> const &trial)
	{
		tmp::d_call_each<
			typename function_space_traits<TestSpace>::field_type_vector_t,
			typename function_space_traits<TrialSpace>::field_type_vector_t,
			eval_on<Operator, tmp::_1, tmp::_2>,
			result_t &,
			Operator const &,
			function_space_base<TestSpace> const &,
			function_space_base<TrialSpace> const &
		>(result, op, test, trial);

		return result;
	}
};

#endif // ASSEMBLY_HPP_INCLUDED

