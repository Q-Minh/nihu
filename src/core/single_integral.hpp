/**
* \file single_integral.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class single_integral
*/

#ifndef SINGLE_INTEGRAL_HPP_INCLUDED
#define SINGLE_INTEGRAL_HPP_INCLUDED

#include "global_definitions.hpp"
#include "../util/plain_type.hpp"
#include "../util/product_type.hpp"
#include "gaussian_quadrature.hpp"
#include "field_type_accelerator.hpp"
#include "formalism.hpp"

template <class TestField, class TrialField>
struct single_integral_traits
{
	/** \brief result type of the single integral residual */
	typedef typename plain_type<
		typename product_type<
			typename TestField::nset_t::shape_t,
			Eigen::Transpose<typename TrialField::nset_t::shape_t >
		>::type
	>::type result_t;
};


/**
* \brief single integral over an element
* \tparam TestField type of the test field
* \tparam TrialField type of the trial field
*/
template <class TestField, class TrialField, class Formalism = typename get_formalism<TestField, TrialField>::type>
class single_integral_impl
{
public:
	/** \brief template parameter as nested type */
	typedef TestField test_field_t;
	/** \brief template parameter as nested type */
	typedef TrialField trial_field_t;

	typedef single_integral_traits<test_field_t, trial_field_t> traits_t;
	typedef typename traits_t::result_t result_t;

	/** \brief L-set of the elem */
	typedef typename test_field_t::elem_t::lset_t lset_t;

	/** \brief the quadrature family */
	typedef gauss_family_tag quadrature_family_t;

	/** \brief N-set of the test field */
	typedef typename test_field_t::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename trial_field_t::nset_t trial_nset_t;

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

		typedef store<field_type_accelerator_pool<
				test_field_t, quadrature_family_t, acceleration::hard, GLOBAL_MAX_ORDER
		> > test_store_t;

		typedef store<field_type_accelerator_pool<
				trial_field_t, quadrature_family_t, acceleration::hard, GLOBAL_MAX_ORDER
		> > trial_store_t;

		auto acc = create_dual_field_type_accelerator(
			test_store_t::m_data[degree],
			trial_store_t::m_data[degree],
			iteration::diagonal());

		for (auto it = acc.begin(); it != acc.end(); ++it)
		{
			auto jac = field.get_elem().get_normal(it.get_first()->get_xi()).norm();
			result += it.get_first()->get_N() *
				(it.get_first()->get_w()*jac) *
				it.get_second()->get_N().transpose();
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
class single_integral_impl<TestField, TrialField, formalism::collocational>
{
public:
	/** \brief N-set of the test field */
	typedef typename TestField::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename TrialField::nset_t trial_nset_t;

	typedef single_integral_traits<TestField, TrialField> traits_t;
	typedef typename traits_t::result_t result_t;

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
template <class TestField, class TrialField, class = void>
class single_integral
{
public:
	typedef single_integral_traits<TestField, TrialField> traits_t;
	typedef typename traits_t::result_t result_t;

	static constexpr result_t eval(
		field_base<TestField> const &,
		field_base<TrialField> const &)
	{
		return result_t::Zero();
	}
};


template <class TestField, class TrialField>
class single_integral<TestField, TrialField,
	typename std::enable_if<
		std::is_same<
			typename TestField::elem_t,
			typename TrialField::elem_t
		>::type::value
	>::type
> :	public single_integral_impl<TestField, TrialField> {};

#endif // SINGLE_INTEGRAL_HPP

