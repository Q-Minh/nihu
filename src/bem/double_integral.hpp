/**
* \file double_integral.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class double_integral and its specialisations
*/

#ifndef DOUBLE_INTEGRAL_HPP_INCLUDED
#define DOUBLE_INTEGRAL_HPP_INCLUDED

#include "quadrature_pool.hpp"
#include "../util/plain_type.hpp"
#include "kernel.hpp"

/**
* \brief class evaluating double integrals of the weighted residual approach
* \tparam Kernel type of the kernel to integrate
* \tparam Test type of the test field
* \tparam Trial type of the trial field
*/
template <bool isTestDirac, class Kernel, class TestField, class TrialField>
class double_integral_impl
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"Kernel must be derived from kernel_base<Kernel>");

	typedef std::true_type WITH_SINGULARITY_CHECK;
	typedef std::false_type WITHOUT_SINGULARITY_CHECK;

public:
	typedef Kernel kernel_t;		/**< \brief template parameter as nested type */
	typedef TestField test_field_t;		/**< \brief template parameter as nested type */
	typedef TrialField trial_field_t;	/**< \brief template parameter as nested type */

	/** \brief test input type of kernel */
	typedef typename kernel_t::test_input_t test_input_t;	
	/** \brief trial input type of kernel */
	typedef typename kernel_t::trial_input_t trial_input_t;
	/** \brief weighted test input type of kernel */
	typedef typename weighted_kernel_input<test_input_t>::type w_test_input_t;
	/** \brief weighted trial input type of kernel */
	typedef typename weighted_kernel_input<trial_input_t>::type w_trial_input_t;
	/** \brief result type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;

	/** \brief the test elem type */
	typedef typename test_field_t::elem_t test_elem_t;
	/** \brief the trial elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;
	/** \brief the test domain type */
	typedef typename test_elem_t::domain_t test_domain_t;
	/** \brief the trial domain type */
	typedef typename trial_elem_t::domain_t trial_domain_t;

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	/** \brief indicates if kernel is singular and singular accelerators need to be instantiated */
	static bool const is_kernel_singular = kernel_traits<kernel_t>::singularity_order != 0;

	/** \brief the quadrature element */
	typedef typename quadrature_type<
		quadrature_family_t,
		typename trial_field_t::elem_t::domain_t
	>::type::quadrature_elem_t quadrature_elem_t;

	typedef regular_pool_store<test_field_t, quadrature_family_t> test_regular_store_t;
	typedef regular_pool_store<trial_field_t, quadrature_family_t> trial_regular_store_t;

	typedef field_type_accelerator<test_field_t, quadrature_family_t> test_field_type_accelerator_t;
	typedef field_type_accelerator<trial_field_t, quadrature_family_t> trial_field_type_accelerator_t;


	/** \brief N-set of the test field */
	typedef typename test_field_t::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename trial_field_t::nset_t trial_nset_t;

	/** \brief type of test shape function */
	typedef typename test_nset_t::shape_t test_shape_t;
	/** \brief type of trial shape function */
	typedef typename trial_nset_t::shape_t trial_shape_t;

	/** \brief L-set of the test field */
	typedef typename test_field_t::elem_t::lset_t test_lset_t;
	/** \brief L-set of the trial field */
	typedef typename trial_field_t::elem_t::lset_t trial_lset_t;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
		kernel_result_t,
		typename product_type<test_shape_t, Eigen::Transpose<trial_shape_t> >::type
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
	static result_t &eval_on_accelerator(
		result_t &result,
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		test_field_type_accelerator_t const &test_acc,
		field_base<trial_field_t> const &trial_field,
		trial_field_type_accelerator_t const &trial_acc)
	{
		for (auto test_it = test_acc.cbegin(); test_it != test_acc.cend(); ++test_it)
		{
			w_test_input_t test_input(test_field.get_elem(), test_it->get_quadrature_elem().get_xi());
			auto bound = kernel.bind(test_input);
			auto left = (test_it->get_shape() * (test_input.get_jacobian() * test_it->get_quadrature_elem().get_w())).eval();
			for (auto trial_it = trial_acc.cbegin(); trial_it != trial_acc.cend(); ++trial_it)
			{
				w_trial_input_t trial_input(trial_field.get_elem(), trial_it->get_quadrature_elem().get_xi());
				auto right = (trial_it->get_shape().transpose() * (trial_input.get_jacobian() * trial_it->get_quadrature_elem().get_w())).eval();
				result += left * bound.eval(trial_input) * right;
			}
		}

		return result;
	}

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
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field,
		singular_iterator_t begin,
		singular_iterator_t end)
	{
		while (begin != end)
		{
			w_test_input_t test_input(test_field.get_elem(), begin.get_test_quadrature_elem().get_xi());
			w_trial_input_t trial_input(trial_field.get_elem(), begin.get_trial_quadrature_elem().get_xi());

			/** \todo check if lazy evaluation is still faster */
			auto left = (test_field_t::nset_t::eval_shape(begin.get_test_quadrature_elem().get_xi())
				* (test_input.get_jacobian() * begin.get_test_quadrature_elem().get_w())).eval();
			auto right = (trial_field_t::nset_t::eval_shape(begin.get_trial_quadrature_elem().get_xi())
				* (trial_input.get_jacobian() * begin.get_trial_quadrature_elem().get_w())).eval();

			result += left * kernel.eval(test_input, trial_input) * right.transpose();

			++begin;
		}

		return result;
	}

	/** \brief evaluate double integral of a kernel on specific fields without singularity check
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field reference to the test field
	* \param [in] trial_field reference to the trial field
	* \return reference to the stored result
	*/
	static result_t &eval_general(
		WITHOUT_SINGULARITY_CHECK,
		result_t &result,
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field)
	{
		auto &test_ra = test_regular_store_t::m_regular_pool;
		auto &trial_ra = trial_regular_store_t::m_regular_pool;

		// select quadrature
		test_input_t test_center(test_field.get_elem(), test_domain_t::get_center());
		trial_input_t trial_center(trial_field.get_elem(), trial_domain_t::get_center());
		/** \todo why only trial size is important? */
		unsigned degree = kernel.estimate_complexity_interface(test_center, trial_center, trial_field.get_elem().get_linear_size_estimate());

		degree += std::max(test_nset_t::polynomial_order, trial_nset_t::polynomial_order)
			+ std::max(test_lset_t::jacobian_order, trial_lset_t::jacobian_order);

		return eval_on_accelerator(
			result,
			kernel,
			test_field, *(test_ra[degree]),
			trial_field, *(trial_ra[degree]));
	}


	/** \brief evaluate double integral of a kernel on specific fields with singularity check
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field reference to the test field
	* \param [in] trial_field reference to the trial field
	* \return reference to the stored result
	*/
	static result_t &eval_general(
		WITH_SINGULARITY_CHECK,
		result_t &result,
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field)
	{
		typedef accel_store<formalism::general, kernel_t, test_field_t, trial_field_t> acc_store_t;
		auto &sa = acc_store_t::m_singular_accelerator;

		// check singularity
		if (sa.is_singular(test_field, trial_field))
			return eval_singular_on_accelerator(
				result, kernel, test_field, trial_field, sa.begin(), sa.end());

		return eval_general(WITHOUT_SINGULARITY_CHECK(),
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
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field,
		OnSameMesh)
	{
		static bool const sing_check_needed =
			is_kernel_singular && std::is_same<OnSameMesh, std::true_type>::value;

		result_t result;
		result.setZero();	// clear result

		return eval_general(std::integral_constant<bool, sing_check_needed>(),
			result, kernel, test_field, trial_field);
	}
};




template <class Kernel, class TestField, class TrialField>
class double_integral_impl<true, Kernel, TestField, TrialField>
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"Kernel must be derived from kernel_base<Kernel>");

	typedef std::true_type WITH_SINGULARITY_CHECK;
	typedef std::false_type WITHOUT_SINGULARITY_CHECK;

public:
	typedef Kernel kernel_t;		/**< \brief template parameter as nested type */
	typedef TestField test_field_t;		/**< \brief template parameter as nested type */
	typedef TrialField trial_field_t;	/**< \brief template parameter as nested type */

	/** \brief test input type of kernel */
	typedef typename kernel_t::test_input_t test_input_t;	
	/** \brief trial input type of kernel */
	typedef typename kernel_t::trial_input_t trial_input_t;
	/** \brief weighted test input type of kernel */
	typedef typename weighted_kernel_input<test_input_t>::type w_test_input_t;
	/** \brief weighted trial input type of kernel */
	typedef typename weighted_kernel_input<trial_input_t>::type w_trial_input_t;
	/** \brief result type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;

	/** \brief the test elem type */
	typedef typename test_field_t::elem_t test_elem_t;
	/** \brief the trial elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;
	/** \brief the test domain type */
	typedef typename test_elem_t::domain_t test_domain_t;
	/** \brief the trial domain type */
	typedef typename trial_elem_t::domain_t trial_domain_t;

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	/** \brief indicates if kernel is singular and singular accelerators need to be instantiated */
	static bool const is_kernel_singular = kernel_traits<kernel_t>::singularity_order != 0;

	/** \brief the quadrature element */
	typedef typename quadrature_type<
		quadrature_family_t,
		typename trial_field_t::elem_t::domain_t
	>::type::quadrature_elem_t quadrature_elem_t;

	typedef regular_pool_store<trial_field_t, quadrature_family_t> trial_regular_store_t;

	typedef field_type_accelerator<test_field_t, quadrature_family_t> test_field_type_accelerator_t;
	typedef field_type_accelerator<trial_field_t, quadrature_family_t> trial_field_type_accelerator_t;


	/** \brief N-set of the test field */
	typedef typename test_field_t::nset_t test_nset_t;
	/** \brief N-set of the trial field */
	typedef typename trial_field_t::nset_t trial_nset_t;

	/** \brief type of test shape function */
	typedef typename test_nset_t::shape_t test_shape_t;
	/** \brief type of trial shape function */
	typedef typename trial_nset_t::shape_t trial_shape_t;

	/** \brief L-set of the test field */
	typedef typename test_field_t::elem_t::lset_t test_lset_t;
	/** \brief L-set of the trial field */
	typedef typename trial_field_t::elem_t::lset_t trial_lset_t;

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
		kernel_result_t,
		typename product_type<test_shape_t, Eigen::Transpose<trial_shape_t> >::type
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
	static result_t &eval_collocational_on_accelerator(
		result_t &result,
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field,
		trial_field_type_accelerator_t const &trial_acc)
	{
		int row = 0;
		for (auto test_it = test_nset_t::corner_begin(); test_it != test_nset_t::corner_end(); ++test_it)
		{
			test_input_t collocational_point(test_field.get_elem(), *test_it);
			auto bound = kernel.bind(collocational_point);
			for (auto trial_it = trial_acc.cbegin(); trial_it != trial_acc.cend(); ++trial_it)
			{
				w_trial_input_t trial_input(trial_field.get_elem(), trial_it->get_quadrature_elem().get_xi());
				result.row(row) += bound.eval(trial_input) *
					(trial_it->get_shape() * (trial_input.get_jacobian() * trial_it->get_quadrature_elem().get_w()));
			}
			++row;
		}

		return result;
	}


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
	static result_t &eval_collocational_singular_on_accelerator(
		result_t &result,
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field,
		singular_accelerator_t const &sa)
	{
		for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
		{
			quadrature_elem_t qe(test_nset_t::corner_at(idx));
			test_input_t collocational_point(test_field.get_elem(), qe.get_xi());
			auto bound = kernel.bind(collocational_point);
			auto const &quad = sa.get_trial_quadrature(idx);

			for (auto quad_it = quad.begin(); quad_it != quad.end(); ++quad_it)
			{
				w_trial_input_t trial_input(trial_field.get_elem(), quad_it->get_xi());

				result.row(idx) += bound.eval(trial_input) *
					(trial_input.get_jacobian() * quad_it->get_w() *
					trial_field_t::nset_t::eval_shape(quad_it->get_xi()));
			}
		}

		return result;
	}


	/** \brief evaluate single integral of a kernel on specific fields without singularity check
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field reference to the test field
	* \param [in] trial_field reference to the trial field
	* \return reference to the stored result
	*/
	static result_t &eval_collocational(
		WITHOUT_SINGULARITY_CHECK,
		result_t &result,
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field)
	{
		typedef regular_pool_store<trial_field_t, quadrature_family_t> regular_trial_store_t;
		auto &trial_ra = trial_regular_store_t::m_regular_pool;

		quadrature_elem_t qe(test_field_t::elem_t::domain_t::get_center());
		test_input_t test_center(test_field.get_elem(), qe.get_xi());
		trial_input_t trial_center(trial_field.get_elem(), trial_domain_t::get_center());

		/** \todo why only trial size is important? */
		unsigned degree = kernel.estimate_complexity_interface(test_center, trial_center, trial_field.get_elem().get_linear_size_estimate());
		degree += trial_nset_t::polynomial_order + trial_lset_t::jacobian_order;

		return eval_collocational_on_accelerator(
			result,  kernel, test_field, trial_field, *(trial_ra[degree]));
	}


	/** \brief evaluate single integral of a kernel on specific fields with singularity check
	* \param [out] result reference to the integration result matrix
	* \param [in] kernel the kernel to integrate
	* \param [in] test_field reference to the test field
	* \param [in] trial_field reference to the trial field
	* \return reference to the stored result
	*/
	static result_t &eval_collocational(
		WITH_SINGULARITY_CHECK,
		result_t &result,
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field)
	{
		typedef accel_store<formalism::collocational, kernel_t, test_field_t, trial_field_t> acc_store_t;
		typename acc_store_t::singular_accelerator_t &sa = acc_store_t::m_singular_accelerator;

		// check singularity
		if (sa.is_singular(test_field, trial_field))
			return eval_collocational_singular_on_accelerator(
			result, kernel, test_field, trial_field, sa);
		else
			return eval_collocational(WITHOUT_SINGULARITY_CHECK(),
			result, kernel, test_field, trial_field);
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
		kernel_t const &kernel,
		field_base<test_field_t> const &test_field,
		field_base<trial_field_t> const &trial_field,
		OnSameMesh)
	{
		static bool const sing_check_needed =
			is_kernel_singular && std::is_same<OnSameMesh, std::true_type>::value;

		result_t result;
		result.setZero();	// clear result

		return eval_collocational(
			std::integral_constant<bool, sing_check_needed>(),
			result, kernel, test_field, trial_field);
	}
};

template <class Kernel, class TestField, class TrialField>
class double_integral :
	public double_integral_impl<field_traits<TestField>::is_dirac, Kernel, TestField, TrialField>
{
};


#endif

