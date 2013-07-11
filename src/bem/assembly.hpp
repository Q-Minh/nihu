#ifndef ASSEMBLY_HPP_INCLUDED
#define ASSEMBLY_HPP_INCLUDED

#include "../util/matrix_block.hpp"
#include "integral_operator.hpp"

template <class TestSpace, class TrialSpace>
class assembly
{
private:
	/** \brief evaluate weighted residual on homogeneous function spaces
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam TestField the test field type over which integration is performed
	 * \tparam TrialField the trial field type over which integration is performed
	 */
	template <bool is_empty, class Operator, class TestField, class TrialField>
	struct eval_on_impl
	{ struct type {
		/** \brief evaluate weighted residual on homogeneous function spaces
		 * \tparam result_t the type of the result matrix
		 * \param [out] result reference to the result matrix the result block is inserted to
		 * \param [in] kernel the kernel to integrate
		 * \param [in] test_space the test function space
		 * \param [in] trial_space the trial function space
		 */
		template <class result_t>
		void operator() (result_t &result,
			Operator const &op,
			function_space_base<TestSpace> const &test_space,
			function_space_base<TrialSpace> const &trial_space)
		{
			for (auto test_it = test_space.derived().template field_begin<TestField>();
				test_it != test_space.derived().template field_end<TestField>(); ++test_it)
			{
				for (auto trial_it = trial_space.derived().template field_begin<TrialField>();
					trial_it != trial_space.derived().template field_end<TrialField>(); ++trial_it)
				{
					block(result, test_it->get_dofs(), trial_it->get_dofs())
						+= op.eval_on_fields(*test_it, *trial_it);
				}
			}
		}
	};};


	template <class Operator, class TestField, class TrialField>
	struct eval_on_impl<true, Operator, TestField, TrialField>
	{ struct type {
		template <class result_t>
		void operator() (result_t &, Operator const &, TestSpace const &, TrialSpace const &)
		{
		}
	};};

	/** \brief evaluate weighted residual on homogeneous function spaces
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam TestField the test field type over which integration is performed
	 * \tparam TrialField the trial field type over which integration is performed
	 */
	template <class Operator, class TestField, class TrialField>
	struct eval_on : eval_on_impl<
		integral_operator_traits<Operator>::is_local && !std::is_same<typename TestField::elem_t, typename TrialField::elem_t>::value,
		Operator, TestField, TrialField
	> {};


public:
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
			TestSpace const &,
			TrialSpace const &
		>(result, op, test.derived(), trial.derived());

		return result;
	}
};


#endif // ASSEMBLY_HPP_INCLUDED

