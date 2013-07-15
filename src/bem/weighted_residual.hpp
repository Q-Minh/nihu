/**
* \file weighted_residual.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class ::weighted_residual
*/
#ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED
#define WEIGHTED_RESIDUAL_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "function_space.hpp"
#include "projection.hpp"


// forward declaration
template <class LeftDerived, class RightDerived>
class wr_sum;

/** \brief base class of all weighted residual expressions
* \tparam Derived the CRTP derived class
*/
template <class Derived>
class wr_base
{
public:
	NIHU_CRTP_HELPERS

	/** \brief evaluate weighted residual into matrix
	* \tparam result_t the result matrix type
	* \param [out] result the result matrix
	* \return reference to the result
	*/
	template <class result_t>
	result_t &eval(result_t &result) const
	{
		return derived().eval(result);
	}
	
	/** \brief construct product of two weighted residuals
	* \tparam OtherDerived the other weighted residual class type
	* \param [in] otherDerived the other weighted residual object
	* \return a wr_sum object proxy
	*/
	template <class OtherDerived>
	wr_sum<Derived, OtherDerived>
		operator+(wr_base<OtherDerived> const &otherDerived) const
	{
		return wr_sum<Derived, OtherDerived>(*this, otherDerived);
	}
};


/** \brief the sum of two weighted residual expressions
* \tparam LeftDerived the left hand side class
* \tparam RightDerived the right hand side class
*/
template <class LeftDerived, class RightDerived>
class wr_sum :
	public wr_base<wr_sum<LeftDerived, RightDerived> >
{
public:
	/** \brief constructor
	* \param [in] left the left hand side wr reference
	* \param [in] right the right hand side wr reference
	*/
	wr_sum(wr_base<LeftDerived> const &left, wr_base<RightDerived> const &right) :
		m_lhs(left.derived()), m_rhs(right.derived())
	{
	}

	/** \brief evaluate wr sum into result
	* \tparam result_t the result matrix type
	* \param [out] result the result matrix
	* \return reference to the result
	*/
	template <class result_t>
	result_t &eval(result_t &result) const
	{
		m_lhs.eval(result);
		m_rhs.eval(result);
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
	/** \brief constructor from test reference and projection reference
	* \param [in] test reference to the test space
	* \param [in] proj reference to the projection
	*/
	weighted_residual(
		function_space_base<TestSpace> const &test,
		projection_base<Projection> const &proj) :
		m_test(test.derived()), m_proj(proj.derived())
	{
	}

	/** \brief evaluate wr sum into result
	* \tparam result_t the result matrix type
	* \param [out] result the result matrix
	* \return reference to the result
	*/
	template <class result_t>
	result_t &eval(result_t &result) const
	{
		m_proj.test_on_into(m_test, result);
		return result;
	}

private:
	/** \brief the stored test function space reference */
	TestSpace const &m_test;
	/** \brief the stored projection reference */
	/** \todo dangerous if projection is created on the fly */
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

