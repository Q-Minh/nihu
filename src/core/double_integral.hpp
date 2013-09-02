/**
* \file double_integral.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
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
template <class Kernel, class TestField, class TrialField, class = typename get_formalism<TestField, TrialField>::type>
class double_integral
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
	/** \brief indicates if kernel is singular and singular accelerators need to be instantiated */
	static bool const is_kernel_singular = kernel_traits<Kernel>::singularity_order != 0;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
			typename kernel_traits<Kernel>::result_t,
			typename product_type<
				typename TestField::nset_t::shape_t,
				Eigen::Transpose<typename TrialField::nset_t::shape_t>
			>::type
		>::type
	>::type result_t;

protected:
	/** \brief evaluate regular double integral with selected accelerators
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field the test field to integrate on
	* \param [in] test_acc field type accelerator of the test field
	* \param [in] trial_field the trial field to integrate on
	* \param [in] trial_acc field type accelerator of the trial field
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

			auto left = (TestField::nset_t::eval_shape(it.get_first()->get_xi())
				* (test_input.get_jacobian() * it.get_first()->get_w())).eval();
			auto right = (TrialField::nset_t::eval_shape(it.get_second()->get_xi())
				* (trial_input.get_jacobian() * it.get_second()->get_w())).eval();

			result += left * kernel(test_input, trial_input) * right.transpose();
		}

		return result;
	}

public:
	/** \brief evaluate double singular integral with selected singular accelerator
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \param [in] begin begin iteartor of the singular quadrature
	* \param [in] end end iterator of the singular quadrature
	* \return reference to the integration result
	*/
	template <class singular_iterator_t>
	static result_t &eval_singular_on_accelerator(
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

			auto left = (TestField::nset_t::eval_shape(begin.get_test_quadrature_elem().get_xi())
				* (test_input.get_jacobian() * begin.get_test_quadrature_elem().get_w())).eval();
			auto right = (TrialField::nset_t::eval_shape(begin.get_trial_quadrature_elem().get_xi())
				* (trial_input.get_jacobian() * begin.get_trial_quadrature_elem().get_w())).eval();

			result += left * kernel(test_input, trial_input) * right.transpose();

			++begin;
		}

		return result;
	}

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
		case REGULAR:
			break;
		case FACE_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, singularity::face_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case CORNER_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, singularity::corner_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case EDGE_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, singularity::edge_match_type
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
			is_kernel_singular && std::is_same<OnSameMesh, std::true_type>::value;

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
	/** \brief indicates if kernel is singular and singular accelerators need to be instantiated */
	static bool const is_kernel_singular = kernel_traits<Kernel>::singularity_order != 0;

	/** \brief N-set of the test field */
	typedef typename TestField::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename TrialField::nset_t trial_nset_t;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
			typename kernel_traits<Kernel>::result_t,
			typename product_type<
				typename test_nset_t::shape_t,
				Eigen::Transpose<typename trial_nset_t::shape_t>
			>::type
		>::type
	>::type result_t;

protected:
	/** \brief evaluate regular collocational integral with selected trial field accelerator
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \param [in] trial_acc the trial field type accelerator
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
	* \tparam singular_accelerator_t type of the singular quadrature accelerator
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \param [in] sa singular accelerator
	* \return reference to the integration result
	*/
	template <class singular_accelerator_t>
	static result_t &eval_singular_on_accelerator(
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
					TrialField::nset_t::eval_shape(quad_it->get_xi()));
			}
		}

		return result;
	}

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
		case REGULAR:
			break;
		case FACE_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, singularity::face_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case CORNER_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, singularity::corner_match_type
				>::eval(result, kernel, test_field, trial_field, match);
		case EDGE_MATCH:
			return singular_integral_shortcut<
				Kernel, TestField, TrialField, singularity::edge_match_type
				>::eval(result, kernel, test_field, trial_field, match);
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
			is_kernel_singular && std::is_same<OnSameMesh, std::true_type>::value;

		result_t result;
		result.setZero();	// clear result

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
	typedef store<
		singular_accelerator<formalism_t, Kernel, TestField, TrialField>
	> store_t;

	template <class result_t>
	static result_t &eval_impl(
		formalism::general,
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &match)
	{
		return double_integral<Kernel, TestField, TrialField>::eval_singular_on_accelerator(
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
		return double_integral<Kernel, TestField, TrialField>::eval_singular_on_accelerator(
			result, kernel, test_field, trial_field, store_t::m_data);
	}


public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result type
	 * \param [out] result the integral result
	 * \param [in] kernel the kernel instance
	 * \param [in] test_field the test field instance
	 * \param [in] trial_field the trial field instance
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

#endif

