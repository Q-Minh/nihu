/**
* \file double_integral.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class double_integral and its specialisations
*/
#ifndef DOUBLE_INTEGRAL_HPP_INCLUDED
#define DOUBLE_INTEGRAL_HPP_INCLUDED

#include <Eigen/Dense>

#include "kernel.hpp"
#include "quadrature.hpp"
#include "kernel_input.hpp"
#include "accelerator.hpp"


/**
 * \brief class evaluating double integrals of the weighted residual approach
 * \tparam Kernel type of the kernel to integrate
 * \tparam Test type of the test field
 * \tparam Trial type of the trial field
 */
template <class Kernel, class Test, class Trial>
class double_integral
{
public:
	typedef Kernel kernel_t;		/**< \brief template parameter as nested type */

	typedef Test test_field_t;		/**< \brief template parameter as nested type */
	typedef Trial trial_field_t;	/**< \brief template parameter as nested type */

	typedef typename kernel_t::input_t kernel_input_t;		/**< \brief input type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;	/**< \brief result type of kernel */

	typedef gauss_quadrature<typename test_field_t::elem_t::domain_t> test_quadrature_t; 	/**< \brief type of test quadrature */
	typedef gauss_quadrature<typename trial_field_t::elem_t::domain_t> trial_quadrature_t;	/**< \brief type of trial quadrature */

	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;	/**< \brief type of quadrature element */

	typedef accelerator_by<typename test_field_t::elem_t> test_accelerator_t;	/**< \brief type of the accelerator of the test field */
	typedef accelerator_by<typename trial_field_t::elem_t> trial_accelerator_t;	/**< \brief type of the accelerator of the trial field */

	typedef typename test_field_t::nset_t test_nset_t;		/**< \brief type of element's N-set */
	typedef typename trial_field_t::nset_t trial_nset_t;	/**< \brief type of element's N-set */

	typedef typename Eigen::Matrix<
		typename kernel_t::scalar_t, test_field_t::num_dofs, trial_field_t::num_dofs
	> result_t;		/**< \brief integration result type */


	/** \brief evaluate double integral with selected quadratures
	 * \param [in] test_field the test field to integrate on
	 * \param [in] test_quad the test quadrature
	 * \param [in] trial_field the trial field to integrate on
	 * \param [in] trial_quad the trial quadrature
	 * \return reference to the integration result
	 */
	static result_t const &eval_on_quadrature(
		test_field_t const &test_field,
		test_quadrature_t const &test_quad,
		trial_field_t const &trial_field,
		trial_quadrature_t const &trial_quad)
	{
		m_result = result_t();	// clear result

		for(auto test_it = test_quad.begin(); test_it != test_quad.end(); ++test_it)
		{
			kernel_input_t test_input(test_field.get_elem(), *test_it);
			for(auto trial_it = trial_quad.begin(); trial_it != trial_quad.end(); ++trial_it)
			{
				kernel_input_t trial_input(trial_field.get_elem(), *trial_it);

				m_result +=
					test_nset_t::eval_shape(test_it->get_xi()) *
					kernel_t::eval(test_input, trial_input) *
					test_input.get_jacobian() *
					trial_input.get_jacobian() *
					trial_nset_t::eval_shape(trial_it->get_xi()).transpose();
			}
		}

		return m_result;
	}

	/** \brief evaluate double integral on given fields
	 * \param [in] test_field the test field to integrate on
	 * \param [in] trial_field the trial field to integrate on
	 * \return reference to the integration result
	 */
	static result_t const &eval(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		auto const &test_quad_pool = m_test_accelerator.get_quadrature_pool();
		auto const &trial_quad_pool = m_trial_accelerator.get_quadrature_pool();

		// select quadrature
		kernel_input_t test_center(test_field.get_elem(), test_quad_pool[0][0]);
		kernel_input_t trial_center(trial_field.get_elem(), trial_quad_pool[0][0]);
		unsigned degree = kernel_t::estimate_complexity(test_center, trial_center);

		return eval_on_quadrature(test_field, test_quad_pool[degree], trial_field, trial_quad_pool[degree]);
	}

protected:
	static result_t m_result; /**< \brief the integral result stored as static variable */

	static const test_accelerator_t m_test_accelerator; /**< \brief accelerator of the test field */
	static const trial_accelerator_t m_trial_accelerator; /**< \brief accelerator of the trial field */
};

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::result_t
	double_integral<Kernel, Test, Trial>::m_result;

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::test_accelerator_t
	const double_integral<Kernel, Test, Trial>::m_test_accelerator;

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::trial_accelerator_t
	const double_integral<Kernel, Test, Trial>::m_trial_accelerator;


/**
 * \brief specialisation of double integral for the collocational approach
 * \tparam Kernel type of the kernel to integrate
 * \tparam Test type of the test field
 * \tparam Trial type of the trial field
 */
template<class Kernel, class Trial, class ElemType, class FieldOption>
class double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>
{
public:
	typedef Kernel kernel_t;		/**< \brief template parameter as nested type */

	typedef field<ElemType, FieldOption, dirac_field> test_field_t;		/**< \brief template parameter as nested type */
	typedef Trial trial_field_t;	/**< \brief template parameter as nested type */

	typedef typename kernel_t::input_t kernel_input_t;		/**< \brief input type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;	/**< \brief result type of kernel */

	typedef gauss_quadrature<typename trial_field_t::elem_t::domain_t> trial_quadrature_t;	/**< \brief type of trial quadrature */
	typedef typename trial_quadrature_t::quadrature_elem_t quadrature_elem_t;	/**< \brief type of quadrature element */

	typedef accelerator_by<typename trial_field_t::elem_t> trial_accelerator_t;	/**< \brief type of the accelerator of the trial field */

	typedef typename test_field_t::nset_t test_nset_t;		/**< \brief type of element's N-set */
	typedef typename trial_field_t::nset_t trial_nset_t;	/**< \brief type of element's N-set */

	typedef typename Eigen::Matrix<
		typename kernel_t::scalar_t, test_field_t::num_dofs, trial_field_t::num_dofs
	> result_t;		/**< \brief integration result type */


	/** \brief evaluate double integral with selected trial field quadrature
	 * \param [in] test_field the test field to integrate on
	 * \param [in] trial_field the trial field to integrate on
	 * \param [in] trial_quad the trial quadrature
	 * \return reference to the integration result
	 */
	static result_t const &eval_on_quadrature(
		test_field_t const &test_field,
		trial_field_t const &trial_field,
		trial_quadrature_t const &trial_quad)
	{
		m_result = result_t();	// clear result

		for(auto test_it = test_nset_t::corner_begin(); test_it != test_nset_t::corner_end(); ++test_it)
		{
			quadrature_elem_t qe(*test_it);
			kernel_input_t collocational_point(test_field.get_elem(), qe);
			for(auto trial_it = trial_quad.begin(); trial_it != trial_quad.end(); ++trial_it)
			{
				kernel_input_t trial_input(trial_field.get_elem(), *trial_it);

				m_result.row(test_it - test_nset_t::corner_begin()) +=
					kernel_t::eval(collocational_point, trial_input) *
					trial_input.get_jacobian() *
					trial_nset_t::eval_shape(trial_it->get_xi()).transpose();
			}
		}

		return m_result;
	}

	/** \brief evaluate double integral on given fields
	 * \param [in] test_field the test field to integrate on
	 * \param [in] trial_field the trial field to integrate on
	 * \return reference to the integration result
	 */
	static result_t const &eval(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		auto const &trial_quad_pool = m_trial_accelerator.get_quadrature_pool();

		// select quadrature
		quadrature_elem_t qe(test_field_t::elem_t::domain_t::get_center());
		kernel_input_t test_center(test_field.get_elem(), qe);
		kernel_input_t trial_center(trial_field.get_elem(), trial_quad_pool[0][0]);
		unsigned degree = kernel_t::estimate_complexity(test_center, trial_center);
		return eval_on_quadrature(test_field, trial_field, trial_quad_pool[degree]);
	}

protected:
	static result_t m_result; /**< \brief the integral result stored as static variable */

	static const trial_accelerator_t m_trial_accelerator; /**< \brief accelerator of the trial field */
};

template<class Kernel, class Trial, class ElemType, class FieldOption>
typename double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::result_t
	double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::m_result;

template<class Kernel, class Trial, class ElemType, class FieldOption>
typename double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::trial_accelerator_t
	const double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::m_trial_accelerator;


#endif

