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
 * \file complexity_estimator.hpp
 * \ingroup kernel
 * \brief Estimate kernel complexity between two elements
 */
#ifndef COMPLEXITY_ESTIMATOR_HPP_INCLUDED
#define COMPLEXITY_ESTIMATOR_HPP_INCLUDED

#include "../tmp/integer.hpp"
#include "field.hpp"

/** \brief clas to estimate kernel complexity between two fields
 * \tparam TestField the test field
 * \tparam TrialField the trial field
 * \tparam KernelEstimator the kernel's complexity estimator
 */
template <class TestField, class TrialField, class KernelEstimator>
class complexity_estimator
{
public:
	/** \brief the test field type */
	typedef TestField test_field_t;
	/** \brief the trial field type */
	typedef TrialField trial_field_t;

	enum {
		/** \brief the test field complexity */
		test_field_complexity =
			shape_set_traits::polynomial_order<typename test_field_t::nset_t>::value +
			shape_set_traits::jacobian_order<typename test_field_t::elem_t::lset_t>::value,
		/** \brief the trial field complexity */
		trial_field_complexity =
			shape_set_traits::polynomial_order<typename trial_field_t::nset_t>::value +
			shape_set_traits::jacobian_order<typename trial_field_t::elem_t::lset_t>::value
	};

	/** \brie the total field complexity */
	static unsigned const total_field_complexity = tmp::max_<
		tmp::integer<unsigned, test_field_complexity>,
		tmp::integer<unsigned, trial_field_complexity>
	>::type::value;

	/** \brief compute total complexity
	 * \param [in] test_field the test field instance
	 * \param [in] trial_field the trial field instance
	 * \return the kernel complexity
	 */
	static unsigned eval(
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field
	)
	{
		return total_field_complexity + KernelEstimator::eval(test_field, trial_field);
	}
};


/**
 * \brief merge at least two complexity estimators
 * \tparam Estim1 the first estimator
 * \tparam Estims the further complexity estimators
 */
template <class Estim1, class...Estims>
struct merge_kernel_complexity_estimators :
	merge_kernel_complexity_estimators<
		Estim1,
		typename merge_kernel_complexity_estimators<Estims...>::type
	>
{
};


/**
 * \brief merge two complexity estimators (the general case)
 * \tparam Estim1 the first estimator
 * \tparam Estim2 the second estimator
 */
template <class Estim1, class Estim2>
struct merge_kernel_complexity_estimators<Estim1, Estim2>
{
	struct type
	{
		/** \brief compute complexity of merged estimators
		 * \tparam test_field_t the test field type
		 * \tparam trial_field_t the trial field type
		 * \param [in] test_field the test field instance
		 * \param [in] trial_field the trial field instance
		 * \return the kernel complexity
		 */
		template <class test_field_t, class trial_field_t>
		static unsigned eval(
			field_base<test_field_t> const &test_field,
			field_base<trial_field_t> const &trial_field
		)
		{
			return std::max(
				Estim1::eval(test_field, trial_field),
				Estim2::eval(test_field, trial_field)
			);
		}
	};
};

#endif

