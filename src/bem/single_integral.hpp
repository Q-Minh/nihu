/**
* \file single_integral.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class single_integral
*/

#ifndef SINGLE_INTEGRAL_HPP_INCLUDED
#define SINGLE_INTEGRAL_HPP_INCLUDED

#include "../util/plain_type.hpp"
#include "../util/product_type.hpp"
#include "gaussian_quadrature.hpp"
#include "quadrature_pool.hpp"

/**
* \brief single integral over an element
* \tparam TestField type of the test field
* \tparam TrialField type of the trial field
*/
template <bool isTestDirac, class TestField, class TrialField>
class single_integral_impl
{
public:
	/** \brief template parameter as nested type */
	typedef TestField test_field_t;
	/** \brief template parameter as nested type */
	typedef TrialField trial_field_t;

	/** \brief the elem type */
	typedef typename test_field_t::elem_t elem_t;
	/** \brief the test domain type */
	typedef typename elem_t::domain_t domain_t;
	/** \brief L-set of the elem */
	typedef typename elem_t::lset_t lset_t;

	/** \brief the quadrature family */
	typedef gauss_family_tag quadrature_family_t;

	/** \brief the stored regular quadratures for the test field */
	typedef regular_pool_store<test_field_t, quadrature_family_t> test_regular_store_t;
	/** \brief the stored regular quadratures for the trial field */
	typedef regular_pool_store<trial_field_t, quadrature_family_t> trial_regular_store_t;

	/** \brief N-set of the test field */
	typedef typename test_field_t::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename trial_field_t::nset_t trial_nset_t;

	/** \brief type of test shape function */
	typedef typename test_nset_t::shape_t test_shape_t;
	/** \brief type of trial shape function */
	typedef typename trial_nset_t::shape_t trial_shape_t;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
			test_shape_t,
			typename plain_type<
				Eigen::Transpose<trial_shape_t>
			>::type
		>::type
	>::type result_t;
	
	/** \brief the polynomial degree of the product of the n-sets and the jacobian */
	static unsigned const degree
		= test_nset_t::polynomial_order
		+ trial_nset_t::polynomial_order
		+ lset_t::jacobian_order;

public:
	/** \brief evaluate single integral on a given field
	* \param [in] field the field to integrate on
	* \return the integration result by value
	*/
	static result_t eval(
		field_base<test_field_t> const &field,
		field_base<trial_field_t> const &)
	{
		result_t result;
		result.setZero();

		auto &test_ra = test_regular_store_t::m_regular_pool;
		auto &trial_ra = trial_regular_store_t::m_regular_pool;
		auto const &test_acc = *(test_ra[degree]);
		auto const &trial_acc = *(trial_ra[degree]);

		elem_t const &elem = field.get_elem();

		auto trial_it = trial_acc.cbegin();
		for (auto test_it = test_acc.cbegin(); test_it != test_acc.cend(); ++test_it, ++trial_it)
		{
			auto xi = test_it->get_quadrature_elem().get_xi();
			auto jac = elem.get_normal(xi).norm();
			auto w = test_it->get_quadrature_elem().get_w();
			result += test_it->get_shape() * (w*jac) * trial_it->get_shape().transpose();
		}

		return result;
	}
};


/**
* \brief specialisation for the collocational case
* \tparam TestField type of the test field
* \tparam TrialField type of the trial field
*/
template <class TestField, class TrialField>
class single_integral_impl<true, TestField, TrialField>
{
public:
	/** \brief N-set of the test field */
	typedef typename TestField::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename TrialField::nset_t trial_nset_t;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
			typename  test_nset_t::shape_t,
			typename plain_type<
				Eigen::Transpose<typename trial_nset_t::shape_t>
			>::type
		>::type
	>::type result_t;
	
public:
	/** \brief evaluate collocational integral on a given field
	* \return integration result by value
	*/
	static result_t eval(
		field_base<TestField> const &,
		field_base<TrialField> const &)
	{
		result_t result;
		result.setZero();
		for (unsigned row = 0; row < test_nset_t::num_nodes; ++row)
			result.row(row) += trial_nset_t::eval_shape(test_nset_t::corner_at(row));
		return result;
	}
};


/**
* \brief single integral over a field
* \tparam TestField type of the test field
* \tparam TrialField type of the trial field
*/
template <class TestField, class TrialField>
class single_integral :
	public single_integral_impl<field_traits<TestField>::is_dirac, TestField, TrialField>
{
};


#endif // SINGLE_INTEGRAL_HPP

