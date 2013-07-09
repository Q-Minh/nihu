#ifndef PROJECTION_HPP_INCLUDED
#define PROJECTION_HPP_INCLUDED

template <class LDerived, class RDerived>
class projection_sum;

template <class Derived>
class projection_base
{
public:
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

	template <class TestSpace, class Result>
	void test_on_into(
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		derived().test_on_into(test_space, result);
	}

	template <class RDerived>
	projection_sum<Derived, RDerived>
		operator+(projection_base<RDerived> const &rhs) const
	{
		return projection_sum<Derived, RDerived>(*this, rhs);
	}
};


template <class LDerived, class RDerived>
class projection_sum : public projection_base<projection_sum<LDerived, RDerived> >
{
public:
	projection_sum(projection_base<LDerived> const &left, projection_base<RDerived> const &right) :
		m_lhs(left.derived()), m_rhs(right.derived())
	{
	}

	template <class TestSpace, class Result>
	void test_on_into(
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		m_lhs.test_on_into(test_space, result);
		m_rhs.test_on_into(test_space, result);
	}
	
private:
	LDerived const m_lhs;
	RDerived const m_rhs;
};



template <class Operator, class TrialSpace>
class projection : public projection_base<projection<Operator, TrialSpace> >
{
public:
	projection(integral_operator_base<Operator> const &op,
		function_space_base<TrialSpace> const &trial) :
		m_op(op.derived()), m_trial_space(trial.derived())
	{
	}

	template <class TestSpace, class Result>
	void test_on_into(
		function_space_base<TestSpace> const &test_space,
		Result &result) const
	{
		std::cout << "Tesztelek, geexci" << std::endl;
	}
	

private:
	Operator m_op;
	TrialSpace const &m_trial_space;
};

template <class Operator, class TrialSpace>
projection<Operator, TrialSpace> operator*(integral_operator_base<Operator> const &op, function_space_base<TrialSpace> const &trial)
{
	return projection<Operator, TrialSpace>(op, trial);
}

#endif

