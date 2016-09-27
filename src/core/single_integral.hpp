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
* \file single_integral.hpp
* \ingroup intop
* \brief declaration of class single_integral
*/

#ifndef SINGLE_INTEGRAL_HPP_INCLUDED
#define SINGLE_INTEGRAL_HPP_INCLUDED

#include "global_definitions.hpp"
#include "../util/store_pattern.hpp"
#include "../util/plain_type.hpp"
#include "../util/product_type.hpp"
#include "../util/block_product.hpp"
#include "gaussian_quadrature.hpp"
#include "field_type_accelerator.hpp"
#include "formalism.hpp"

/** \brief traits class of ::single_integral
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
struct single_integral_traits
{
	/** \brief result type of the single integral residual */
	typedef typename block_product_result_type<
		typename TestField::nset_t::shape_t,
		Eigen::Matrix<double, TestField::quantity_dimension, TrialField::quantity_dimension>,
		typename TrialField::nset_t::shape_t
	>::type result_t;
};


/**
* \brief single integral over an element for the general case
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

	/** \brief the traits class */
	typedef single_integral_traits<test_field_t, trial_field_t> traits_t;
	/** \brief the result matrix type */
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
				test_field_t, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > test_store_t;

		typedef store<field_type_accelerator_pool<
				trial_field_t, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > trial_store_t;

		auto acc = create_dual_field_type_accelerator(
			test_store_t::get_data()[degree],
			trial_store_t::get_data()[degree],
			iteration::diagonal());
			
		// identity matrix instance for block product integration
		enum {Rows = TestField::quantity_dimension, Cols = TrialField::quantity_dimension};
		auto mat = Eigen::Matrix<double, Rows, Cols>::Identity();

		for (auto it = acc.begin(); it != acc.end(); ++it)
		{
			auto jac = field.get_elem().get_normal(it.get_first()->get_xi()).norm();
			result += block_product(it.get_first()->get_N(),
				mat * it.get_first()->get_w()*jac,
				it.get_second()->get_N());
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

	/** \brief the traits class */
	typedef single_integral_traits<TestField, TrialField> traits_t;
	/** \brief the result matrix type of the single integral */
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
		
		enum {Rows = TestField::quantity_dimension, Cols = TrialField::quantity_dimension};
		auto mat = Eigen::Matrix<double, Rows, Cols>::Identity();
		
		for (unsigned row = 0; row < test_nset_t::num_nodes; ++row)
		{
			auto N = trial_nset_t::template eval_shape<0>(test_nset_t::corner_at(row));
			for (unsigned col = 0; col < trial_nset_t::num_nodes; ++col)
				result.template block<Rows, Cols>(row*Rows, col*Cols) = mat * N(col);
		}
		
		return result;
	}
};


/**
 * \brief single integral for different element types
 * \tparam TestField type of the test field
 * \tparam TrialField type of the trial field
 */
template <class TestField, class TrialField, class = void>
class single_integral
{
public:
	/** \brief the result matrix type */
	typedef typename single_integral_traits<TestField, TrialField>::result_t result_t;

	/** \brief specialisation of single_integral::eval for the empty case
	 * \return the empty result matrix
	 */
	static constexpr result_t eval(
		field_base<TestField> const &,
		field_base<TrialField> const &)
	{
		return result_t::Zero();
	}
};


/**
 * \brief single integral for matching element types
 * \tparam TestField type of the test field
 * \tparam TrialField type of the trial field
 */
template <class TestField, class TrialField>
class single_integral<TestField, TrialField,
	typename std::enable_if<
		std::is_same<
			typename TestField::elem_t,
			typename TrialField::elem_t
		>::value
	>::type
> :	public single_integral_impl<TestField, TrialField> {};

#endif // SINGLE_INTEGRAL_HPP

