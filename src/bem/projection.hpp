/**
* \file projection.hpp
* \ingroup intop
* \brief declaration of projection classes
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/

#ifndef PROJECTION_HPP_INCLUDED
#define PROJECTION_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "function_space.hpp"
#include "integral_operator.hpp"
#include "assembly.hpp"

/** \brief CRTP base class of all projection expressions
* \tparam Derived CRTP derived class
* \details a projection is an integral operator multiplied by a function space.
*/
template <class Derived>
class projection_base
{
public:
	NIHU_CRTP_HELPERS

	/** \brief test a projection with a test function space
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


/** \brief metafunction determining if argument is projection expression
 * \tparam Proj the function space to test
 */
template <class Proj>
struct is_projection : std::is_base_of<
	projection_base<typename std::decay<Proj>::type>,
	typename std::decay<Proj>::type
>{};


/** \brief proxy class representing a sum of two projections
* \tparam LDerived left hand side projection expression
* \tparam RDerived right hand side projection expression
*/
template <class LDerived, class RDerived>
class projection_sum :
	public projection_base<projection_sum<LDerived, RDerived> >
{
public:
	/** \brief constructor from two projection expressions
	* \param [in] left the left hand side projection expression
	* \param [in] right the right hand side projection expression
	*/
	projection_sum(LDerived &&left, RDerived &&right) :
		m_lhs(std::forward<LDerived>(left)),
		m_rhs(std::forward<RDerived>(right))
	{
	}

	/** \brief test a projection with a test function space
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
	/** \brief the left hand side projection expression */
	LDerived m_lhs;
	/** \brief the right hand side projection expression */
	RDerived m_rhs;
};


/** \brief factory operator for the sum of two projections */
template <class Left, class Right>
projection_sum<
	typename std::enable_if<is_projection<Left>::value, Left>::type,
	typename std::enable_if<is_projection<Right>::value, Right>::type
>
	operator+(Left &&left, Right &&right)
{
	return projection_sum<Left, Right>(
		std::forward<Left>(left),
		std::forward<Right>(right));
}


/** \brief integral operator applied on a function space
* \tparam Operator the integral operator
* \tparam TrialSpace the right hand side function space
*/
template <class Operator, class TrialSpace>
class projection
	: public projection_base<projection<Operator, TrialSpace> >
{
public:
	/** \brief constructor from an operator and a function space
	* \param [in] op the operator
	* \param [in] trial the trial function space
	*/
	projection(integral_operator_base<Operator> const &op, TrialSpace && trial) :
		m_op(op.derived()),
		m_trial_space(std::forward<TrialSpace>(trial))
	{
	}


private:
	template <class TestSpace, class Result>
	void test_on_into(
		std::true_type,
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		if (&test_space.get_mesh() == &m_trial_space.get_mesh())
			assembly<TestSpace, typename std::decay<TrialSpace>::type, std::true_type>::eval_into(
				result, m_op, test_space, m_trial_space);
		else
			assembly<TestSpace, typename std::decay<TrialSpace>::type, std::false_type>::eval_into(
				result, m_op, test_space, m_trial_space);
	}

	template <class TestSpace, class Result>
	void test_on_into(
		std::false_type,
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		assembly<TestSpace, typename std::decay<TrialSpace>::type, std::false_type>::eval_into(
			result, m_op, test_space, m_trial_space);
	}

public:
	/** \brief test a projection with a test function space
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
		test_on_into(
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

