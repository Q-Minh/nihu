// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file assembly.hpp
 * \ingroup assembly
 * \brief assemble WR matrices from field wr submatrices
 */
#ifndef ASSEMBLY_HPP_INCLUDED
#define ASSEMBLY_HPP_INCLUDED

#include "../util/matrix_block.hpp"
#include "integral_operator.hpp"
#include "function_space.hpp"

/** \brief assemble result matrix from field wr submatrices
* \tparam TestSpace the test function space
* \tparam TrialSpace the trial function space
*/
template <class TestSpace, class TrialSpace, class OnSameMesh>
class assembly
{
private:
	/** \brief implementation of ::eval_on
	 * \tparam Operator the integral operator type that is applied
	 * \tparam TestField the test field type over which integration is performed
	 * \tparam TrialField the trial field type over which integration is performed
	 * \tparam isTrivial indicates if the operation is void or not
	 */
	template <class Operator, class TestField, class TrialField, bool isTrivial =
		integral_operator_traits<Operator>::is_local &&
		!std::is_same<typename TestField::elem_t, typename TrialField::elem_t>::value
		>
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
			integral_operator_base<Operator> const &op,
			function_space_base<TestSpace> const &test_space,
			function_space_base<TrialSpace> const &trial_space)
		{
			if (!integral_operator_traits<Operator>::is_local)
			{
				for (auto test_it = test_space.template field_begin<TestField>();
					test_it != test_space.template field_end<TestField>(); ++test_it)
					for (auto trial_it = trial_space.template field_begin<TrialField>();
						trial_it != trial_space.template field_end<TrialField>(); ++trial_it)
						block(result, test_it->get_dofs(), trial_it->get_dofs())
							+= op.eval_on_fields(*test_it, *trial_it, OnSameMesh());
			}
			else	// local operator
			{
				// we can be sure that the element types are the same and element id's are monotonic
				auto test_it = test_space.template field_begin<TestField>();
				auto test_end = test_space.template field_end<TestField>();
				auto trial_it = trial_space.template field_begin<TrialField>();
				auto trial_end = trial_space.template field_end<TrialField>();

				for (; test_it != test_end; ++test_it)
				{
					auto test_id = test_it->get_elem().get_id();
					for (; trial_it != trial_end && trial_it->get_elem().get_id() < test_id; ++trial_it);
					if (trial_it->get_elem().get_id() == test_id)
						block(result, test_it->get_dofs(), trial_it->get_dofs())
							+= op.eval_on_fields(*test_it, *trial_it, OnSameMesh());
				}
			}
		}
	};};


	/** \brief trivial specialisation of ::eval_on_impl for the trivial case */
	template <class Operator, class TestField, class TrialField>
	struct eval_on_impl<Operator, TestField, TrialField, true>
	{ struct type {
		template <class result_t>
		void operator() (
			result_t &,
			integral_operator_base<Operator> const &,
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
	struct eval_on : eval_on_impl<Operator, TestField, TrialField>
	{
	};


public:
	/** \brief evaluate weighted residual into result matrix
	* \tparam result_t the result matrix type
	* \tparam Operator the operator type
	* \param [out] result the result matrix reference
	* \param [in] op reference to the integral operator
	* \param [in] test_space reference to the test function space
	* \param [in] trial_space reference to the trial function space
	*/
	template <class result_t, class Operator>
	static result_t &eval_into(
		result_t &result,
		integral_operator_base<Operator> const &op,
		function_space_base<TestSpace> const &test_space,
		function_space_base<TrialSpace> const &trial_space)
	{
		if(!OnSameMesh::value && integral_operator_traits<Operator>::is_local)
			return result;

		tmp::d_call_each<
			typename function_space_traits<TestSpace>::field_type_vector_t,
			typename function_space_traits<TrialSpace>::field_type_vector_t,
			eval_on<Operator, tmp::_1, tmp::_2>,
			result_t &,
			integral_operator_base<Operator> const &,
			function_space_base<TestSpace> const &,
			function_space_base<TrialSpace> const &
		>(result, op, test_space, trial_space);

		return result;
	}
};

#endif // ASSEMBLY_HPP_INCLUDED

