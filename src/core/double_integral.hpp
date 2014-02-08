// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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
* \file double_integral.hpp
* \ingroup intop
* \brief declaration of class double_integral and its specialisations
*/

#ifndef DOUBLE_INTEGRAL_HPP_INCLUDED
#define DOUBLE_INTEGRAL_HPP_INCLUDED

#include "global_definitions.hpp"
#include "../util/product_type.hpp"
#include "../util/plain_type.hpp"
#include "../util/brick.hpp"
#include "../library/location_normal.hpp"
#include "kernel.hpp"
#include "complexity_estimator.hpp"
#include "element_match.hpp"
#include "field_type_accelerator.hpp"
#include "singular_accelerator.hpp"
#include "formalism.hpp"

// forward declaration
template <class Kernel, class TestField, class TrialField, class Singularity, class Enable = void>
class singular_integral_shortcut;

/**
* \brief class evaluating double integrals of the weighted residual approach
* \tparam Kernel type of the kernel to integrate
* \tparam TestField type of the test field
* \tparam TrialField type of the trial field
*/
template <
	class Kernel, class TestField, class TrialField,
	class Formalism = typename get_formalism<TestField, TrialField>::type
>
class double_integral;


/**
 * \brief specialisation of ::double_integral for the general formalism
 * \tparam Kernel type of the kernel to integrate
 * \tparam TestField type of the test field
 * \tparam TrialField type of the trial field
 */
template <class Kernel, class TestField, class TrialField>
class double_integral<Kernel, TestField, TrialField, formalism::general>
{
	typedef std::true_type WITH_SINGULARITY_CHECK;
	typedef std::false_type WITHOUT_SINGULARITY_CHECK;

public:
	/** \brief test input type of kernel */
	typedef typename kernel_traits<Kernel>::test_input_t test_input_t;
	/** \brief trial input type of kernel */
	typedef typename kernel_traits<Kernel>::trial_input_t trial_input_t;
	/** \brief weighted test input type of kernel */
	typedef typename merge<
		test_input_t,
		typename build<normal_jacobian<typename test_input_t::space_t> >::type
	>::type w_test_input_t;
	/** \brief weighted trial input type of kernel */
	typedef typename merge<
		trial_input_t,
		typename build<normal_jacobian<typename trial_input_t::space_t> >::type
	>::type w_trial_input_t;

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<Kernel>::quadrature_family_t quadrature_family_t;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
			typename kernel_traits<Kernel>::result_t,
			typename plain_type<
				typename product_type<
					typename TestField::nset_t::shape_t,
					Eigen::Transpose<typename TrialField::nset_t::shape_t>
				>::type
			>::type
		>::type
	>::type result_t;

protected:
	/** \brief evaluate regular double integral with selected accelerators
	 * \param [out] result reference to the integration result matrix
	 * \param [in] kernel the kernel to integrate
	 * \param [in] test_field the test field to integrate on
	 * \param [in] trial_field the trial field to integrate on
	 * \param [in] it begin iterator of the accelerator
	 * \param [in] end end iterator of the accelerator
	 * \return reference to the integration result
	 */
	template <class dual_iterator_t>
	static result_t &eval_on_accelerator(
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		dual_iterator_t it,
		dual_iterator_t end
		)
	{
		for (; it != end; ++it)
		{
			w_test_input_t test_input(test_field.get_elem(), it.get_first()->get_xi());
			w_trial_input_t trial_input(trial_field.get_elem(), it.get_second()->get_xi());

			auto left = (TestField::nset_t::template eval_shape<0>(it.get_first()->get_xi())
				* (test_input.get_jacobian() * it.get_first()->get_w())).eval();
			auto right = (TrialField::nset_t::template eval_shape<0>(it.get_second()->get_xi())
				* (trial_input.get_jacobian() * it.get_second()->get_w())).eval();

			result += left * kernel(test_input, trial_input) * right.transpose();
		}

		return result;
	}

public:
	/** \brief evaluate double singular integral with selected singular accelerator
	 * \tparam singular_accelerator_t the singular accelerator class
	 * \tparam dummy a dummy type needed to keep specialisation within the class body
	 */
	template <class singular_accelerator_t, class dummy>
	struct eval_singular_on_accelerator
	{
		/** \brief evaluate double singular integral with selected singular accelerator
		 * \param [out] result reference to the integration result matrix
		 * \param [in] kernel the kernel to integrate
		 * \param [in] test_field the test field to integrate on
		 * \param [in] trial_field the trial field to integrate on
		 * \param [in] begin begin iterator of the singular quadrature
		 * \param [in] end end iterator of the singular quadrature
		 * \return reference to the integration result
		 */
		template <class singular_iterator_t>
		static result_t &eval(
			result_t &result,
			kernel_base<Kernel> const &kernel,
			field_base<TestField> const &test_field,
			field_base<TrialField> const &trial_field,
			singular_iterator_t begin,
			singular_iterator_t end)
		{
			while (begin != end)
			{
				w_test_input_t test_input(test_field.get_elem(), begin.get_test_quadrature_elem().get_xi());
				w_trial_input_t trial_input(trial_field.get_elem(), begin.get_trial_quadrature_elem().get_xi());

				auto left = (TestField::nset_t::template eval_shape<0>(begin.get_test_quadrature_elem().get_xi())
					* (test_input.get_jacobian() * begin.get_test_quadrature_elem().get_w())).eval();
				auto right = (TrialField::nset_t::template eval_shape<0>(begin.get_trial_quadrature_elem().get_xi())
					* (trial_input.get_jacobian() * begin.get_trial_quadrature_elem().get_w())).eval();

				result += left * kernel(test_input, trial_input) * right.transpose();

				++begin;
			}

			return result;
		}
	};

	/** \brief specialisation of eval_singular_on_accelerator to the invalid accelerator
	 * \tparam dummy a dummy type needed to keep specialisation within the class body
	 */
	template <class dummy>
	struct eval_singular_on_accelerator<invalid_singular_accelerator, dummy>
	{
		/** \brief evaluate double singular integral with the invalid accelerator
		 * \details this function simply throws an exception
		 */
		template <class singular_iterator_t>
		static result_t &eval(
			result_t &result,
			kernel_base<Kernel> const &,
			field_base<TestField> const &,
			field_base<TrialField> const &,
			singular_iterator_t,
			singular_iterator_t)
		{
			throw std::runtime_error("Invalid general singular quadrature");
			return result;
		}
	};

protected:
	/** \brief evaluate double integral of a kernel on specific fields without singularity check
	 * \param [out] result reference to the integration result matrix
	 * \param [in] kernel the kernel to integrate
	 * \param [in] test_field reference to the test field
	 * \param [in] trial_field reference to the trial field
	 * \return reference to the stored result
	 */
	static result_t &eval(
		WITHOUT_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
		unsigned degree = complexity_estimator<
			TestField, TrialField,
			typename kernel_traits<Kernel>::complexity_estimator_t
		>::eval(test_field, trial_field);

		typedef store<field_type_accelerator_pool<
			TestField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > test_store_t;

		typedef store<field_type_accelerator_pool<
			TrialField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > trial_store_t;

		auto acc = create_dual_field_type_accelerator(
			test_store_t::m_data[degree], trial_store_t::m_data[degree], iteration::diadic());

		return eval_on_accelerator(
			result, kernel, test_field, trial_field, acc.begin(), acc.end());
	}


	/** \brief evaluate double integral of a kernel on specific fields with singularity check
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field reference to the test field
	* \param [in] trial_field reference to the trial field
	* \return reference to the stored result
	*/
	static result_t &eval(
		WITH_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
		auto match = element_match_eval(test_field, trial_field);
		switch (match.get_singularity_type())
		{
		case NO_MATCH:
			break;
		case FACE_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, match::face_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case CORNER_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, match::corner_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case EDGE_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, match::edge_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		}
		return eval(WITHOUT_SINGULARITY_CHECK(), result, kernel, test_field, trial_field);
	}

public:
	/** \brief evaluate double integral on given fields
	 * \param [in] kernel the kernel to integrate
	 * \param [in] test_field the test field to integrate on
	 * \param [in] trial_field the trial field to integrate on
	 * \return the integration result by value
	 */
	template <class OnSameMesh>
	static result_t eval(
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		OnSameMesh)
	{
		static bool const sing_check_needed =
			kernel_traits<Kernel>::is_singular && std::is_same<OnSameMesh, std::true_type>::value;

		result_t result;
		result.setZero();	// clear result

		return eval(std::integral_constant<bool, sing_check_needed>(),
			result, kernel, test_field, trial_field);
	}
};




/**
* \brief specialisation of ::double_integral for the collocational formalism
* \tparam Kernel type of the kernel to integrate
* \tparam TestField type of the test field
* \tparam TrialField type of the trial field
*/
template <class Kernel, class TestField, class TrialField>
class double_integral<Kernel, TestField, TrialField, formalism::collocational>
{
	typedef std::true_type WITH_SINGULARITY_CHECK;
	typedef std::false_type WITHOUT_SINGULARITY_CHECK;

public:
	/** \brief test input type of kernel */
	typedef typename kernel_traits<Kernel>::test_input_t test_input_t;
	/** \brief trial input type of kernel */
	typedef typename kernel_traits<Kernel>::trial_input_t trial_input_t;

	/** \brief weighted test input type */
	typedef test_input_t w_test_input_t;
	/** \brief weighted trial input type */
	typedef typename merge<
		trial_input_t,
		typename build<normal_jacobian<typename trial_input_t::space_t> >::type
	>::type w_trial_input_t;

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<Kernel>::quadrature_family_t quadrature_family_t;

	/** \brief N-set of the test field */
	typedef typename TestField::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename TrialField::nset_t trial_nset_t;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
			typename kernel_traits<Kernel>::result_t,
			typename plain_type<	// this plain type is needed by clang
				typename product_type<
					typename test_nset_t::shape_t,
					Eigen::Transpose<typename trial_nset_t::shape_t>
				>::type
			>::type
		>::type
	>::type result_t;

protected:
	/** \brief evaluate regular collocational integral with selected trial field accelerator
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \param [in] it the begin iterator of the accelerator
	* \param [in] end the end iterator of the accelerator
	* \return reference to the integration result
	*/
	template <class dual_iterator_t>
	static result_t &eval_on_accelerator(
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		dual_iterator_t it,
		dual_iterator_t end)
	{
		for (; it != end; ++it)
		{
			w_test_input_t test_input(test_field.get_elem(), it.get_first()->get_xi());
			w_trial_input_t trial_input(trial_field.get_elem(), it.get_second()->get_xi());

			auto left = (it.get_first()->get_N() * (it.get_first()->get_w())).eval();
			auto right = (it.get_second()->get_N() * (trial_input.get_jacobian() * it.get_second()->get_w())).eval();

			result += left * kernel(test_input, trial_input) * right.transpose();
		}

		return result;
	}

public:
	/** \brief evaluate collocational singular integral with selected singular accelerator
	 * \tparam singular_accelerator_t the singular accelerator type
	 * \tparam dummy dummy argument to keep explicit specialisation within class body
	 */
	template <class singular_accelerator_t, class dummy>
	struct eval_singular_on_accelerator
	{
		/** \brief evaluate collocational singular integral with selected singular accelerator
		* \tparam singular_accelerator_t type of the singular quadrature accelerator
		* \param [out] result reference to the integration result matrix
		* \param [in] kernel the kernel to integrate
		* \param [in] test_field the test field to integrate on
		* \param [in] trial_field the trial field to integrate on
		* \param [in] sa singular accelerator
		* \return reference to the integration result
		*/
		static result_t & eval(
			result_t &result,
			kernel_base<Kernel> const &kernel,
			field_base<TestField> const &test_field,
			field_base<TrialField> const &trial_field,
			singular_accelerator_t const &sa)
		{
			for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
			{
				test_input_t collocational_point(test_field.get_elem(), test_nset_t::corner_at(idx));
				auto bound = kernel.bind(collocational_point);
				auto const &quad = sa.get_trial_quadrature(idx);

				for (auto quad_it = quad.begin(); quad_it != quad.end(); ++quad_it)
				{
					w_trial_input_t trial_input(trial_field.get_elem(), quad_it->get_xi());

					result.row(idx) += bound(trial_input) *
						(trial_input.get_jacobian() * quad_it->get_w() *
						TrialField::nset_t::template eval_shape<0>(quad_it->get_xi()));
				}
			}

			return result;
		}
	};

	/** \brief evaluate collocational singular integral with the invalid accelerator
	 * \tparam dummy dummy argument to keep explicit specialisation within class body
	 */
	template <class dummy>
	struct eval_singular_on_accelerator<invalid_singular_accelerator, dummy>
	{
		/** \brief evaluate collocational singular integral with selected singular accelerator
		* \tparam singular_accelerator_t type of the singular quadrature accelerator
		* \param [out] result reference to the integration result matrix
		* \return reference to the integration result
		* \details throws and exception is called, this case does not exist.
		*/
		static result_t & eval(
			result_t &result,
			kernel_base<Kernel> const &,
			field_base<TestField> const &,
			field_base<TrialField> const &,
			invalid_singular_accelerator const &)
		{
			throw std::runtime_error("Invalid quadrature");
			return result;
		}
	};

protected:
	/** \brief evaluate single integral of a kernel on specific fields without singularity check
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field reference to the test field
	* \param [in] trial_field reference to the trial field
	* \return reference to the stored result
	*/
	static result_t &eval(
		WITHOUT_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
		unsigned degree = complexity_estimator<
			TestField, TrialField,
			typename kernel_traits<Kernel>::complexity_estimator_t
		>::eval(test_field, trial_field);

		if (degree > GLOBAL_MAX_ORDER)
			throw std::out_of_range("Too high quadrature degree selected for collocational integration");

		typedef store<field_type_accelerator_pool<
			TestField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > test_store_t;

		typedef store<field_type_accelerator_pool<
			TrialField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > trial_store_t;

		auto acc = create_dual_field_type_accelerator(
			test_store_t::m_data[degree], trial_store_t::m_data[degree], iteration::diadic());

		return eval_on_accelerator(
			result, kernel, test_field, trial_field, acc.begin(), acc.end());
	}


	/** \brief evaluate single integral of a kernel on specific fields with singularity check
	 * \param [out] result reference to the integration result matrix
	 * \param [in] kernel the kernel to integrate
	 * \param [in] test_field reference to the test field
	 * \param [in] trial_field reference to the trial field
	 * \return reference to the stored result
	 */
	static result_t &eval(
		WITH_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
		auto match = element_match_eval(test_field, trial_field);
		switch (match.get_singularity_type())
		{
		case NO_MATCH:
			break;
		case FACE_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, match::face_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case CORNER_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, match::corner_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case EDGE_MATCH:	// just for the warning, this case is impossible
			;
		}
		return eval(WITHOUT_SINGULARITY_CHECK(), result, kernel, test_field, trial_field);
	}

public:
	/** \brief evaluate collocational integral on given fields
	* \tparam singularity_check_needed indicates if surface system is solved or not
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \return the integration result by value
	*/
	template <class OnSameMesh>
	static result_t eval(
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		OnSameMesh)
	{
		static bool const sing_check_needed =
			kernel_traits<Kernel>::is_singular && std::is_same<OnSameMesh, std::true_type>::value;

		result_t result;
		result.setZero();

		return eval(
			std::integral_constant<bool, sing_check_needed>(),
			result, kernel, test_field, trial_field);
	}
};


/** \brief a shortcut for the user to define customised singular integral methods
 * \tparam Kernel the kernel class
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 * \tparam Singularity the singularity type
 * \tparam Enable additional argument for std::enable_if
 */
template <class Kernel, class TestField, class TrialField, class Singularity, class Enable>
class singular_integral_shortcut
{
	typedef typename get_formalism<TestField, TrialField>::type formalism_t;
	typedef typename select_singular_accelerator<Kernel, TestField, TrialField>::type singular_accelerator_t;
	typedef store<singular_accelerator_t> store_t;

	template <class result_t>
	static result_t &eval_impl(
		formalism::general,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &match)
	{
		return double_integral<Kernel, TestField, TrialField>::template eval_singular_on_accelerator<singular_accelerator_t, void>::eval(
			result, kernel, test_field, trial_field, store_t::m_data.begin(match), store_t::m_data.end(match));
	}

	template <class result_t>
	static result_t &eval_impl(
		formalism::collocational,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		return double_integral<Kernel, TestField, TrialField>::template eval_singular_on_accelerator<singular_accelerator_t, void>::eval(
			result, kernel, test_field, trial_field, store_t::m_data);
	}

public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result type
	 * \param [out] result the integral result
	 * \param [in] kernel the kernel instance
	 * \param [in] test_field the test field instance
	 * \param [in] trial_field the trial field instance
	 * \param [in] match the element match information
	 * \return reference to the integral result
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &match)
	{
		return eval_impl(formalism_t(), result, kernel, test_field, trial_field, match);
	}
};


/** \brief trivial overload for singular_integral_shortcut for the collocation with constant test field
 * \todo this specialisation should not be implemented, try to avoid referencing
 * \tparam Kernel the kernel class
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 * \tparam Singularity the singularity type
 * \tparam Enable additional argument for std::enable_if
 */
template <class Kernel, class TestField, class TrialField, class Singularity>
class singular_integral_shortcut<Kernel, TestField, TrialField, Singularity,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TestField::nset_t, constant_shape_set<typename TestField::lset_t::domain_t> >::value &&
		std::is_same<Singularity, match::corner_match_type>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result type
	 * \param [out] result the integral result
	 * \return reference to the integral result
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<Kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};

#endif

