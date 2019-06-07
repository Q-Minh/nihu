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
* \file double_integral.hpp
* \ingroup intop
* \brief declaration of class double_integral and its specialisations
*/

#ifndef DOUBLE_INTEGRAL_HPP_INCLUDED
#define DOUBLE_INTEGRAL_HPP_INCLUDED

#include "../util/matrix_traits.hpp"
#include "../util/product_type.hpp"
#include "../util/plain_type.hpp"
#include "../util/brick.hpp"
#include "../util/block_product.hpp"
#include "../library/location_normal.hpp"
#include "../tmp/control.hpp"
#include "../util/store_pattern.hpp"
#include "complexity_estimator.hpp"
#include "element_match.hpp"
#include "field_type_accelerator.hpp"
#include "formalism.hpp"
#include "global_definitions.hpp"
#include "kernel.hpp"
#include "match_types.hpp"
#include "nearly_singular_integral.hpp"
#include "singular_accelerator.hpp"
#include "singular_integral_shortcut.hpp"

namespace NiHu
{


/*
template <class Elem>
struct weighted_brick;

template <class LSet, class Scalar>
struct weighted_brick<volume_element<LSet, Scalar> >
{
	typedef volume_jacobian<typename volume_element<LSet, Scalar>::space_t> type;
};

template <class LSet, class Scalar>
struct weighted_brick<surface_element<LSet, Scalar> >
{
	typedef normal_jacobian<typename surface_element<LSet, Scalar>::space_t> type;
};
*/


/** Helper class for tmp::call_until to select singular shortcuts based on singularity type */
template <class Singularity>
struct singular_shortcut_switch
{
	struct type
	{
		template <class result_t, class Kernel, class TestField, class TrialField>
		bool operator()(
			result_t &result,
			kernel_base<Kernel> const &kernel,
			field_base<TestField> const &test_field,
			field_base<TrialField> const &trial_field,
			element_match const &mtch)
		{
#if NIHU_DEBUGGING
			std::cout << "Match dimension: " << mtch.get_match_dimension() << std::endl;
			std::cout << "Singularity value: " << Singularity::value << std::endl;
#endif

#if NIHU_MEX_DEBUGGING
			static bool printed = false;
			if (!printed)
			{
				mexPrintf("Singular shortcut switch called for mtch dim: %d, sing val: %d\n", mtch.get_match_dimension(), Singularity::value);
				printed = true;
			}
#endif

			// if the parameter singularity is valid, evaluate shortcut
			if (mtch.get_match_dimension() == Singularity::value)
			{
				singular_integral_shortcut<Kernel, TestField, TrialField, Singularity>::eval(
					result, kernel, test_field, trial_field, mtch);
				return true;
			}

			return false;
		}
	};
};



template <class Kernel, class TestField, class TrialField>
struct double_integral_traits
{
	/** \brief result type of the weighted residual */
	typedef typename block_product_result_type<
		typename TestField::nset_t::shape_t,
		typename Kernel::result_t,
		typename TrialField::nset_t::shape_t
	>::type result_t;
};



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
 * \brief specialisation of NiHu::double_integral for the general formalism
 * \tparam Kernel type of the kernel to integrate
 * \tparam TestField type of the test field
 * \tparam TrialField type of the trial field
 */
template <class Kernel, class TestField, class TrialField>
class double_integral<Kernel, TestField, TrialField, formalism::general>
{
	typedef double_integral_traits<Kernel, TestField, TrialField> traits_t;

public:
	typedef std::true_type WITH_SINGULARITY_CHECK;
	typedef std::false_type WITHOUT_SINGULARITY_CHECK;

	/** \brief the test elem type */
	typedef typename TestField::elem_t test_elem_t;
	/** \brief the trial elem type */
	typedef typename TrialField::elem_t trial_elem_t;
	/** \brief test input type of kernel */
	typedef typename kernel_traits<Kernel>::test_input_t test_input_t;
	/** \brief trial input type of kernel */
	typedef typename kernel_traits<Kernel>::trial_input_t trial_input_t;
	/** \brief weighted test input type of kernel */
	typedef typename weighted_input<test_input_t, test_elem_t>::type w_test_input_t;
	/** \brief weighted trial input type of kernel */
	typedef typename weighted_input<trial_input_t, trial_elem_t>::type w_trial_input_t;

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<Kernel>::quadrature_family_t quadrature_family_t;

	/** \brief result type of the weighted residual */
	typedef typename traits_t::result_t result_t;

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

			result += block_product(
				TestField::nset_t::template eval_shape<0>(it.get_first()->get_xi())	*
				(test_input.get_jacobian() * it.get_first()->get_w()),
				kernel(test_input, trial_input),
				TrialField::nset_t::template eval_shape<0>(it.get_second()->get_xi()) *
				(trial_input.get_jacobian() * it.get_second()->get_w())
			);
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
			for (; begin != end; ++begin)
			{
				w_test_input_t test_input(test_field.get_elem(), begin.get_test_quadrature_elem().get_xi());
				w_trial_input_t trial_input(trial_field.get_elem(), begin.get_trial_quadrature_elem().get_xi());

				result += block_product(
					TestField::nset_t::template eval_shape<0>(begin.get_test_quadrature_elem().get_xi()) *
					(test_input.get_jacobian() * begin.get_test_quadrature_elem().get_w()),
					kernel(test_input, trial_input),
					TrialField::nset_t::template eval_shape<0>(begin.get_trial_quadrature_elem().get_xi()) *
					(trial_input.get_jacobian() * begin.get_trial_quadrature_elem().get_w())
				);
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
#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("double_integral::eval without singularity check called on elements %d <- %d\n",
				test_field.get_dofs()(0), trial_field.get_dofs()(0));
			printed = true;
		}
#endif
		// store type of the regular test quadratures
		typedef store<field_type_accelerator_pool<
			TestField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > test_store_t;

		// store type of the regular trial quadratures
		typedef store<field_type_accelerator_pool<
			TrialField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > trial_store_t;

		// determine integration degree based on the complexity estimator
		unsigned degree = complexity_estimator<
			TestField, TrialField,
			typename Kernel::estimator_t
		>::eval(test_field, trial_field);

		auto acc(create_dual_field_type_accelerator(
			test_store_t::get_data()[degree], trial_store_t::get_data()[degree], iteration::diadic()));

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
		// evaluate element match
		auto mtch(element_match_eval(test_field, trial_field));

		// if no match simple regular integral is computed
		if (mtch.get_match_dimension() == -1)
			return eval(WITHOUT_SINGULARITY_CHECK(), result, kernel, test_field, trial_field);

		typedef typename match_type_vector<TestField, TrialField>::type possible_match_types;

		// traverse possible singular integral shortcuts with tmp::call_until
		if (!tmp::call_until<
			possible_match_types,
			singular_shortcut_switch<tmp::_1>,
			result_t &,
			kernel_base<Kernel> const &,
			field_base<TestField> const &,
			field_base<TrialField> const &,
			element_match const &
		>(result, kernel, test_field, trial_field, mtch))
			std::cerr << "UNHANDLED GALERKIN SINGULARITY TYPE: " << mtch.get_match_dimension() << std::endl;

		return result;
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
* \brief specialisation of NiHu::double_integral for the collocational formalism
* \tparam Kernel type of the kernel to integrate
* \tparam TestField type of the test field
* \tparam TrialField type of the trial field
*/
template <class Kernel, class TestField, class TrialField>
class double_integral<Kernel, TestField, TrialField, formalism::collocational>
{

	typedef double_integral_traits<Kernel, TestField, TrialField> traits_t;

public:
	typedef std::true_type WITH_SINGULARITY_CHECK;
	typedef std::false_type WITHOUT_SINGULARITY_CHECK;

	/** \brief type of the trial element */
	typedef typename TrialField::elem_t trial_elem_t;
	/** \brief test input type of kernel */
	typedef typename kernel_traits<Kernel>::test_input_t test_input_t;
	/** \brief trial input type of kernel */
	typedef typename kernel_traits<Kernel>::trial_input_t trial_input_t;
	/** \brief result type of kernel */
	typedef typename kernel_traits<Kernel>::result_t kernel_result_t;

	/** \brief weighted trial input type */
	typedef typename weighted_input<trial_input_t, trial_elem_t>::type w_trial_input_t;

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<Kernel>::quadrature_family_t quadrature_family_t;

	/** \brief N-set of the test field */
	typedef typename TestField::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename TrialField::nset_t trial_nset_t;

	/** \brief result type of the weighted residual */
	typedef typename traits_t::result_t result_t;

	enum {
		/// \brief number of rows of the kernel result
		kernel_rows = num_rows<kernel_result_t>::value,
		/// \brief number of columns of the kernel result
		kernel_cols = num_cols<kernel_result_t>::value
	};

protected:
	/** \brief evaluate regular collocational integral with selected trial field accelerator
	 * \param [out] result reference to the integration result matrix
	 * \param [in] kernel the kernel to integrate
	 * \param [in] test_field defining the collocation points
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
			test_input_t test_input(test_field.get_elem(), it.get_first()->get_xi());
			w_trial_input_t trial_input(trial_field.get_elem(), it.get_second()->get_xi());

			result += block_product(
				it.get_first()->get_N() * (it.get_first()->get_w()),
				kernel(test_input, trial_input),
				it.get_second()->get_N() * (trial_input.get_jacobian() * it.get_second()->get_w())
			);
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
		 * \param [out] result reference to the integration result matrix
		 * \param [in] kernel the kernel to integrate
		 * \param [in] test_field the test field defining the collocation points
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
				test_input_t collocation_point(test_field.get_elem(), test_nset_t::corner_at(idx));

				for (auto const &q : sa.get_trial_quadrature(idx))
				{
					w_trial_input_t trial_input(trial_field.get_elem(), q.get_xi());

					result.template block<kernel_rows, kernel_cols*trial_nset_t::num_nodes>(idx*kernel_rows, 0)
						+= semi_block_product(
							kernel(collocation_point, trial_input),
							trial_nset_t::template eval_shape<0>(q.get_xi()) * trial_input.get_jacobian() * q.get_w()
						);
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
		 * \param [out] result reference to the integration result matrix
		 * \return reference to the integration result
		 * \details throws and exception if called, this case does not exist.
		 */
		static result_t & eval(
			result_t &result,
			kernel_base<Kernel> const &,
			field_base<TestField> const &,
			field_base<TrialField> const &,
			invalid_singular_accelerator const &)
		{
			throw std::runtime_error("Invalid quadrature returned by eval_singular_on_accelerator");
			return result;
		}
	};

protected:
	/** \brief evaluate single integral of a kernel on specific fields without singularity check
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field reference to the test field defining the collocation points
	* \param [in] trial_field reference to the trial field to integrate on
	* \return reference to the stored result
	*/
	static result_t &eval(
		WITHOUT_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
		// evaluate nearly with nearly singular shortcut if needed
		if (nearly_singular_integral<Kernel, TestField, TrialField>::needed(
			kernel, test_field, trial_field))
			return nearly_singular_integral<Kernel, TestField, TrialField>::eval(
				result, kernel, test_field, trial_field);

		// evaluate regular integral
		unsigned degree = complexity_estimator<
			TestField, TrialField,
			typename Kernel::estimator_t
		>::eval(test_field, trial_field);

#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("double_integral::eval without singularity check called on elements %d <- %d\n",
				test_field.get_dofs()(0), trial_field.get_dofs()(0));
			mexPrintf("Regular integral with degree %d\n", degree);
			printed = true;
		}
#endif

		if (degree > GLOBAL_MAX_ORDER)
			throw std::out_of_range("Too high quadrature degree selected for collocational integration");

		typedef store<field_type_accelerator_pool<
			TestField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > test_store_t;

		typedef store<field_type_accelerator_pool<
			TrialField, quadrature_family_t, GLOBAL_ACCELERATION, GLOBAL_MAX_ORDER
		> > trial_store_t;

		auto acc(create_dual_field_type_accelerator(
			test_store_t::get_data()[degree], trial_store_t::get_data()[degree], iteration::diadic()));

		return eval_on_accelerator(
			result, kernel, test_field, trial_field, acc.begin(), acc.end());
	}


	/** \brief evaluate single integral of a kernel on specific fields with singularity check
	 * \param [out] result reference to the integration result matrix
	 * \param [in] kernel the kernel to integrate
	 * \param [in] test_field reference to the test field defining the collocation points
	 * \param [in] trial_field reference to the trial field to integrate on
	 * \return reference to the stored result
	 */
	static result_t &eval(
		WITH_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{

#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("double_integral<Collocation>::eval_with_sing_check called.\n");
			printed = true;
		}
#endif

		auto mtch(element_match_eval(test_field, trial_field));
		if (mtch.get_match_dimension() == -1)
			return eval(WITHOUT_SINGULARITY_CHECK(), result, kernel, test_field, trial_field);

		typedef typename match_type_vector<TestField, TrialField>::type possible_match_types;

#if NIHU_DEBUGGING
		std::cout << "Now checking singular shortcuts" << std::endl;
#endif

		if (!tmp::call_until<
			possible_match_types,
			singular_shortcut_switch<tmp::_1>,
			result_t &,
			kernel_base<Kernel> const &,
			field_base<TestField> const &,
			field_base<TrialField> const &,
			element_match const &
		>(result, kernel, test_field, trial_field, mtch))
		{
			std::cerr << "UNHANDLED COLLOCATIONAL SINGULARITY TYPE: " << mtch.get_match_dimension() << std::endl;
#if NIHU_MEX_DEBUGGING
			mexPrintf("UNHANDLED COLLOCATIONAL SINGULARITY TYPE %d", mtch.get_match_dimension());
#endif

		}
		return result;
	}

public:
	/** \brief evaluate collocational integral on given fields
	 * \tparam OnSameMesh indicates that the two fields are defined over the same mesh and singularity check may be needed
	 * \param [in] kernel the kernel to integrate
	 * \param [in] test_field the test field defining the collocation points
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

#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("double_integral<collocation>::eval called\nSing needed: %d\n", sing_check_needed);
			printed = true;
		}
#endif

		result_t result = result_t::Zero();

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
public:
	typedef void unspecialised;

private:

	typedef typename get_formalism<TestField, TrialField>::type formalism_t;
	typedef typename select_singular_accelerator<
		Kernel, TestField, TrialField, Singularity
	>::type singular_accelerator_t;
	typedef store<singular_accelerator_t> store_t;
	typedef double_integral<Kernel, TestField, TrialField> double_integral_t;

	template <class result_t>
	static result_t &eval_impl(
		formalism::general,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &match)
	{
		return double_integral_t::template eval_singular_on_accelerator<singular_accelerator_t, void>::eval(
			result, kernel, test_field, trial_field,
			store_t::get_data().begin(match),
			store_t::get_data().end(match));
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
#if NIHU_DEBUGGING
		std::cout << "General version of singular_integral_shortcut called" << std::endl;
#endif
#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("General singular integral shortcut for collocation called\n");
			printed = true;
		}
#endif
		return double_integral_t::template eval_singular_on_accelerator<singular_accelerator_t, void>::eval(
			result, kernel, test_field, trial_field, store_t::get_data());
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
#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("General singular integral shortcut called\n");
			printed = true;
		}
#endif
		return eval_impl(formalism_t(), result, kernel, test_field, trial_field, match);
	}
};

}


#endif
