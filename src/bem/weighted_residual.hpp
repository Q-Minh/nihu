/**
* \file weighted_residual.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class ::weighted_residual
*/
#ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED
#define WEIGHTED_RESIDUAL_HPP_INCLUDED

#include "function_space.hpp"
#include "double_integral.hpp"
#include "../util/dual_range.hpp"


template <class LeftDerived, class RightDerived>
class wr_sum;

template <class Derived>
class wr_base
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

	template <class result_t>
	result_t &eval(result_t &result) const
	{
		return derived().eval(result);
	}
	
	template <class OtherDerived>
	wr_sum<Derived, OtherDerived>
		operator+(wr_base<OtherDerived> const &otherDerived) const
	{
		return wr_sum<Derived, OtherDerived>(*this, otherDerived);
	}
};


template <class LeftDerived, class RightDerived>
class wr_sum :
	public wr_base<wr_sum<LeftDerived, RightDerived> >
{
public:
	wr_sum(wr_base<LeftDerived> const &left, wr_base<RightDerived> const &right) :
		m_lhs(left.derived()), m_rhs(right.derived())
	{
	}

	template <class result_t>
	result_t &eval(result_t &result) const
	{
		m_lhs.eval_into(result);
		m_rhs.eval_into(result);
		return result;
	}
	
private:
	LeftDerived const &m_lhs;
	RightDerived const &m_rhs;
};


template <class TestSpace, class Projection>
class weighted_residual :
	public wr_base<weighted_residual<TestSpace, Projection> >
{
public:
	weighted_residual(function_space_base<TestSpace> const &test, projection_base<Projection> const &proj) :
		m_test(test.derived()), m_proj(proj.derived())
	{
	}

	template <class result_t>
	result_t &eval(result_t &result) const
	{
		m_proj.test_on_into(m_test, result);
		return result;
	}

private:
	/** \todo dangerous if projection is created on the fly */
	TestSpace const &m_test;
	Projection const &m_proj;
};

/** \brief factory operator to create a weighted residual from a test space and a projection
* \tparam Test the test space type
* \tparam Proj the projection type
* \param [in] test the test space reference
* \param [in] proj the projection reference
* \return the weighted residual
*/
template <class Test, class Proj>
weighted_residual<Test, Proj>
	operator *(function_space_base<Test> const &test, projection_base<Proj> const &proj)
{
	return weighted_residual<Test, Proj>(test, proj);
}

#endif // ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED

