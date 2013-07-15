/**
* \file projection.hpp
* \ingroup intop
* \brief declaration of projection classes
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/

#ifndef PROJECTION_HPP_INCLUDED
#define PROJECTION_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "assembly.hpp"

// forward declaration
template <class LDerived, class RDerived>
class projection_sum;

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

	/** \brief add a projection to an other (factory operator)
	* \tparam RDerived the right hand side projection type
	* \param [in] rhs the right hand side projection
	* \return the projection sum proxy */
	template <class RDerived>
	projection_sum<Derived, RDerived>
		operator+(projection_base<RDerived> const &rhs) const
	{
		return projection_sum<Derived, RDerived>(*this, rhs);
	}
};


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
	projection_sum(projection_base<LDerived> const &left, projection_base<RDerived> const &right) :
		m_lhs(left.derived()), m_rhs(right.derived())
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
	LDerived const m_lhs;
	/** \brief the right hand side projection expression */
	RDerived const m_rhs;
};



/** \brief proxy class representing the product of an integral operator and a function space
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
	projection(integral_operator_base<Operator> const &op,
		function_space_base<TrialSpace> const &trial) :
		m_op(op.derived()), m_trial_space(trial.derived())
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
		if (&test_space.get_mesh() == &m_trial_space.get_mesh())
			assembly<TestSpace, TrialSpace, std::true_type>::eval_into(
				result, m_op, test_space, m_trial_space);
		else
			assembly<TestSpace, TrialSpace, std::false_type>::eval_into(
				result, m_op, test_space, m_trial_space);
	}
	

private:
	/** \brief the operator stored by value */
	Operator m_op;
	/** \brief the function space stored by reference */
	TrialSpace const &m_trial_space;
};

#endif

