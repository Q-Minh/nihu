/**
* \file weighted_residual.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class ::weighted_residual
*/
#ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED
#define WEIGHTED_RESIDUAL_HPP_INCLUDED

#include <type_traits>

#include "../tmp/bool.hpp"
#include "../util/crtp_base.hpp"
#include "result_matrix.hpp"
#include "function_space.hpp"
#include "integral_transform.hpp"


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
};


/** \brief a weighted residual \f$\left<u,\mathcal{K}v\right>\f$
 * \tparam TestSpace the test function space
 * \tparam Projection the itnegral_transform
 */
template <class TestSpace, class Projection>
class weighted_residual :
	public wr_base<weighted_residual<TestSpace, Projection> >
{
public:
	/** \brief constructor from test reference and integral_transform reference
	* \param [in] test reference to the test space
	* \param [in] proj reference to the integral_transform
	*/
	weighted_residual(
		TestSpace &&test,
		Projection &&proj) :
		m_test(std::forward<TestSpace>(test)),
		m_proj(std::forward<Projection>(proj))
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
	TestSpace m_test;
	/** \brief the stored integral_transform reference */
	/** \todo dangerous if prointegral_transformjection is created on the fly */
	Projection m_proj;
};


/** \brief metafunction determining if artument is weighted_residual expression
 * \tparam Wr the expression to examine
 */
template <class Wr>
struct is_weighted_residual : std::is_base_of<
	wr_base<typename std::decay<Wr>::type>,
	typename std::decay<Wr>::type
>{};


/** \brief the sum of two weighted residual expressions
* \tparam LeftDerived the left hand side class
* \tparam RightDerived the right hand side class
*/
template <class Left, class Right>
class wr_sum :
	public wr_base<wr_sum<Left, Right> >
{
public:
	/** \brief constructor
	* \param [in] left the left hand side wr reference
	* \param [in] right the right hand side wr reference
	*/
	wr_sum(Left &&left, Right &&right) :
		m_lhs(std::forward<Left>(left)), m_rhs(std::forward<Right>(right))
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
	Left m_lhs;
	Right m_rhs;
};


/** \brief factory operator for sum of wr-s
 * \tparam LeftWr left hand side wr
 * \tparam RightWr right hand side wr
 * \param [in] lhs the left hand side
 * \param [in] rhs the right hand side
 */
template <class LeftWr, class RightWr>
wr_sum<
	typename std::enable_if<is_weighted_residual<LeftWr>::value, LeftWr>::type,
	typename std::enable_if<is_weighted_residual<RightWr>::value, RightWr>::type
>
	operator+(LeftWr &&lhs, RightWr &&rhs)
{
	return wr_sum<LeftWr, RightWr>(
		std::forward<LeftWr>(lhs),
		std::forward<RightWr>(rhs));
}



/** \brief factory operator to create a wr from a test space and an integral_transform
* \tparam Test the test space type
* \tparam Proj the integral_transform type
* \param [in] test the test space reference
* \param [in] proj the integral_transform reference
* \return the weighted residual
*/
template <class Test, class Proj>
weighted_residual<
	typename std::enable_if<is_function_space<Test>::value, Test>::type,
	typename std::enable_if<is_integral_transform<Proj>::value, Proj>::type
>
	operator *(Test &&test, Proj &&proj)
{
	return weighted_residual<Test, Proj>(
		std::forward<Test>(test),
		std::forward<Proj>(proj));
}



/** \brief operator to evaluate a weighted residual into a result matrix
 * \tparam WR the weighted residual
 * \tparam Res the result matrix
 * \param [in] wr the weighted residual instance
 * \param [in, out] res the result matrix instance
 * \return (reference to) the result matrix
 */
template <class WR, class Res>
typename std::enable_if<
	is_weighted_residual<WR>::value && is_result_matrix<Res>::value,
	Res &
>::type
operator << (Res &res, WR &&wr)
{
	wr.eval(res);
	return res;
}

/** \brief operator to evaluate a weighted residual into a result matrix
 * \tparam WR the weighted residual
 * \tparam Res the result matrix
 * \param [in] wr the weighted residual instance
 * \param [in, out] res the result matrix instance
 * \return (reference to) the result matrix
 */
template <class WR, class CRes>
typename std::enable_if<
	is_weighted_residual<WR>::value && is_couple<CRes>::value,
	CRes
>::type
operator << (CRes &&res, WR &&wr)
{
	wr.eval(res);
	return res;
}

#endif // ifndef WEIGHTED_RESIDUAL_HPP_INCLUDED

