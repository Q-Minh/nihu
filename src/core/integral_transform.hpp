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
* \file integral_transform.hpp
* \ingroup intop
* \brief declaration of integral_transform classes
*/

#ifndef PROJECTION_HPP_INCLUDED
#define PROJECTION_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "function_space.hpp"
#include "integral_operator.hpp"
#include "assembly.hpp"

/** \brief CRTP base class of all integral_transform expressions
* \tparam Derived CRTP derived class
* \details a integral_transform is an integral operator multiplied by a function space.
*/
template <class Derived>
class integral_transform_base
{
public:
	NIHU_CRTP_HELPERS

	/** \brief test a integral_transform with a test function space
	* \tparam TestSpace type of the test space
	* \tparam Result type of the result matrix
	* \param [in] test_space the test space to test with
	* \param [in, out] result the result matrix to write the result into
	*/
	template <class TestSpace, class Result>
	void test_on_into(
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		derived().test_on_into(test_space, result);
	}
};


/** \brief metafunction determining if argument is integral_transform expression
 * \tparam Proj the function space to test
 */
template <class Proj>
struct is_integral_transform : std::is_base_of<
	integral_transform_base<typename std::decay<Proj>::type>,
	typename std::decay<Proj>::type
>{};


/** \brief proxy class representing a sum of two integral_transforms
* \tparam LDerived left hand side integral_transform expression
* \tparam RDerived right hand side integral_transform expression
*/
template <class LDerived, class RDerived>
class integral_transform_sum :
	public integral_transform_base<integral_transform_sum<LDerived, RDerived> >
{
public:
	/** \brief constructor from two integral_transform expressions
	* \param [in] left the left hand side integral_transform expression
	* \param [in] right the right hand side integral_transform expression
	*/
	integral_transform_sum(LDerived &&left, RDerived &&right) :
		m_lhs(std::forward<LDerived>(left)),
		m_rhs(std::forward<RDerived>(right))
	{
	}

	/** \brief test a integral_transform with a test function space
	* \tparam TestSpace type of the test space
	* \tparam Result type of the result matrix
	* \param [in] test_space the test space to test with
	* \param [in, out] result the result matrix to write the result into
	*/
	template <class TestSpace, class Result>
	void test_on_into(
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		m_lhs.test_on_into(test_space, result);
		m_rhs.test_on_into(test_space, result);
	}

private:
	/** \brief the left hand side integral_transform expression */
	LDerived m_lhs;
	/** \brief the right hand side integral_transform expression */
	RDerived m_rhs;
};


/** \brief factory operator for the sum of two integral_transforms */
template <class Left, class Right>
integral_transform_sum<
	typename std::enable_if<is_integral_transform<Left>::value, Left>::type,
	typename std::enable_if<is_integral_transform<Right>::value, Right>::type
>
	operator+(Left &&left, Right &&right)
{
	return integral_transform_sum<Left, Right>(
		std::forward<Left>(left),
		std::forward<Right>(right));
}


/** \brief Proxy class of an ::integral_operator applied on a ::function_space
* \tparam Operator the integral operator
* \tparam TrialSpace the right hand side function space
*/
template <class Operator, class TrialSpace>
class integral_transform
	: public integral_transform_base<integral_transform<Operator, TrialSpace> >
{
public:
	/** \brief constructor from an operator and a function space
	* \param [in] op the operator
	* \param [in] trial the trial function space
	*/
	integral_transform(integral_operator_base<Operator> const &op, TrialSpace && trial) :
		m_op(op.derived()),
		m_trial_space(std::forward<TrialSpace>(trial))
	{
	}


private:
	template <class TestSpace, class Result>
	void test_on_into_impl(
		std::true_type,
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		if (&test_space.get_mesh() == &m_trial_space.get_mesh())    // singular possible
			assembly<TestSpace, typename std::decay<TrialSpace>::type, std::true_type>::eval_into(
				result, m_op, test_space, m_trial_space);
		else    // singular impossible
			assembly<TestSpace, typename std::decay<TrialSpace>::type, std::false_type>::eval_into(
				result, m_op, test_space, m_trial_space);
	}

	template <class TestSpace, class Result>
	void test_on_into_impl(
		std::false_type,
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		assembly<TestSpace, typename std::decay<TrialSpace>::type, std::false_type>::eval_into(
			result, m_op, test_space, m_trial_space);
	}

public:
	/** \brief test an integral_transform with a test function space
	* \tparam TestSpace type of the test space
	* \tparam Result type of the result matrix
	* \param [in] test_space the test space to test with
	* \param [in, out] result the result matrix to write the result into
	*/
	template <class TestSpace, class Result>
	void test_on_into(
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		test_on_into_impl(
			typename std::is_same<
				typename std::decay<TrialSpace>::type::mesh_t,
				typename std::decay<TestSpace>::type::mesh_t
			>::type(),
			test_space, result);
	}


private:
	/** \brief the operator stored by value */
	Operator m_op;
	/** \brief the function space stored by reference */
	TrialSpace m_trial_space;
};

#endif

