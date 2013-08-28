/**
* \file double_integral.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class double_integral and its specialisations
*/

#ifndef DOUBLE_INTEGRAL_HPP_INCLUDED
#define DOUBLE_INTEGRAL_HPP_INCLUDED

#include "../util/product_type.hpp"
#include "../util/plain_type.hpp"
#include "../util/brick.hpp"
#include "../library/location_normal.hpp"
#include "quadrature_pool.hpp"
#include "kernel.hpp"
#include "complexity_estimator.hpp"
#include "element_match.hpp"

// forward declaration
template <class Kernel, class TestField, class TrialField, class Enable = void>
class singular_integral_shortcut;

/**
* \brief class evaluating double integrals of the weighted residual approach
* \tparam Kernel type of the kernel to integrate
* \tparam Test type of the test field
* \tparam Trial type of the trial field
*/
template <class Kernel, class TestField, class TrialField>
class double_integral
{
public:
	/** \brief template parameter as nested type */
	typedef Kernel kernel_t;
	/** \brief template parameter as nested type */
	typedef TestField test_field_t;
	/** \brief template parameter as nested type */
	typedef TrialField trial_field_t;

	friend class singular_integral_shortcut<Kernel, TestField, TrialField>;

	typedef std::true_type WITH_SINGULARITY_CHECK;
	typedef std::false_type WITHOUT_SINGULARITY_CHECK;

	/** \brief test input type of kernel */
	typedef typename kernel_traits<kernel_t>::test_input_t test_input_t;
	/** \brief trial input type of kernel */
	typedef typename kernel_traits<kernel_t>::trial_input_t trial_input_t;
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

	/** \brief the test elem type */
	typedef typename test_field_t::elem_t test_elem_t;
	/** \brief the trial elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	/** \brief indicates if kernel is singular and singular accelerators need to be instantiated */
	static bool const is_kernel_singular = kernel_traits<kernel_t>::singularity_order != 0;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
			typename kernel_traits<kernel_t>::result_t,
			typename product_type<
				typename test_field_t::nset_t::shape_t,
				Eigen::Transpose<typename trial_field_t::nset_t::shape_t>
			>::type
		>::type
	>::type result_t;

protected:
	template <class iterator_t>
	static result_t &eval_on_accelerator(
		result_t &result,
		kernel_base<kernel_t> const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field,
		iterator_t begin,
		iterator_t end)
	{
		w_test_input_t test_input;
		w_trial_input_t trial_input;

		while (begin != end)
		{
			if (begin.xi_new())
				test_input = test_input(test_field.get_elem(), begin.get_xi());
			if (begin.eta_new())
				trial_input = trial_input(trial_field.get_elem(), begin.get_eta());

			result += (
				begin.get_Nxi() *
				(test_input.get_jacobian() * trial_input.get_jacobian() * begin.get_w()) *
				begin.get_Neta().transpose())
				* kernel(test_input, trial_input);

			++begin;
		}

		return result;
	}


	static result_t &eval(
		WITHOUT_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<kernel_t> const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field)
	{
		typedef regular_pool_store<test_field_t, quadrature_family_t> test_regular_store_t;
		auto &test_ra = test_regular_store_t::m_regular_pool;

		typedef regular_pool_store<trial_field_t, quadrature_family_t> trial_regular_store_t;
		auto &trial_ra = trial_regular_store_t::m_regular_pool;

		unsigned degree = complexity_estimator<
			test_field_t, trial_field_t,
			typename kernel_traits<kernel_t>::complexity_estimator_t
		>::eval(test_field, trial_field);

		auto acc = create_dual_quadrature(
			test_ra[degree].begin(),
			test_ra[degree].end(),
			trial_ra[degree].begin(),
			trial_ra[degree].end(),
			iteration::matrix()
		);

		return eval_on_accelerator(
			result, kernel,
			test_field, trial_field,
			acc.begin(), acc.end());
	}


	static result_t &eval(
		WITH_SINGULARITY_CHECK,
		result_t &result,
		kernel_base<kernel_t> const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field)
	{
		auto match = element_match_eval(test_field, trial_field);
		if (match.get_singularity_type() != REGULAR)
		{
			/** \todo pass singularity type */
			return singular_integral_shortcut<
				kernel_t, test_field_t, trial_field_t>::eval(
				result, kernel, test_field, trial_field);
		}

		return eval(WITHOUT_SINGULARITY_CHECK(),
			result, kernel, test_field, trial_field);
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
		kernel_base<kernel_t> const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field,
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


/** \brief a shortcut for the user to define customised singular integral methods
 * \tparam Formalism collocational or general
 * \tparam Kernel the kernel class
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 * \tparam Enable additional argument for std::enable_if
 */
template <class Kernel, class TestField, class TrialField, class Enable>
class singular_integral_shortcut
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<Kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field)
	{
		typedef singular_store_t<Kernel, TestField, TrialField> sing_store_t;
		auto const &sa = sing_store_t::m_singular_accelerator;

		return double_integral<Kernel, TestField, TrialField>::eval_on_accelerator(
			result, kernel, test_field, trial_field, sa.begin(), sa.end());
	}
};

#endif

