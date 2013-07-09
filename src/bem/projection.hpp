/**
* \file projection.hpp
* \ingroup intop
* \brief declaration of projection classes
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/

#ifndef PROJECTION_HPP_INCLUDED
#define PROJECTION_HPP_INCLUDED

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
	/** \brief CRTP helper function */
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	/** \brief CRTP helper function */
	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
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
	/** \brief reference to the left hand side projection expression */
	LDerived const m_lhs;
	/** \brief reference to the right hand side projection expression */
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
		std::cout << "Tesztelek, geexci" << std::endl;
	}
	

private:
	/** \brief the operator stored by value */
	Operator m_op;
	/** \brief the function space stored by reference */
	TrialSpace const &m_trial_space;
};

/** \brief factory operator to create a projection from an integral operator and a trial space
* \tparam Operator the integral operator type
* \tparam TrialSpace the trial function space type
* \param [in] op the integral operator
* \param [in] trial the trial space reference
* \return projection proxy class
*/
template <class Operator, class TrialSpace>
projection<Operator, TrialSpace>
	operator*(integral_operator_base<Operator> const &op,
	function_space_base<TrialSpace> const &trial)
{
	return projection<Operator, TrialSpace>(op, trial);
}

#endif

