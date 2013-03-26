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
#include "descriptor.hpp"
#include "accelerator.hpp"


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

	typedef accelerator_by<typename test_field_t::elem_t> test_accelerator_t;
	typedef accelerator_by<typename trial_field_t::elem_t> trial_accelerator_t;

	typedef typename test_field_t::nset_t test_nset_t;		/**< \brief type of element's N-set */
	typedef typename trial_field_t::nset_t trial_nset_t;	/**< \brief type of element's N-set */

	typedef typename Eigen::Matrix<
		typename kernel_t::scalar_t, test_field_t::num_dofs, trial_field_t::num_dofs
	> result_t;		/**< \brief integration result type */


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

				m_result += test_nset_t::eval_shape(test_it->get_xi()) *
					kernel_t::eval(test_input, trial_input) *
					test_input.get_jacobian() *
					trial_input.get_jacobian() *
					trial_nset_t::eval_shape(trial_it->get_xi()).transpose();
			}
		}

		return m_result;
	}

	static result_t const &eval(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		auto test_quad_pool = m_test_accelerator.get_quadrature_pool();
		auto trial_quad_pool = m_trial_accelerator.get_quadrature_pool();

		// select quadrature
		kernel_input_t test_center(test_field.get_elem(), test_quad_pool[0][0]);
		kernel_input_t trial_center(trial_field.get_elem(), trial_quad_pool[0][0]);
		unsigned degree = kernel_t::estimate_complexity(test_center, trial_center);

		return eval_on_quadrature(test_field, test_quad_pool[degree], trial_field, trial_quad_pool[degree]);
	}

protected:
	static result_t m_result; /**< \brief the integral result stored as static variable */

	static test_accelerator_t m_test_accelerator; /**< \brief accelerator of the test field */
	static trial_accelerator_t m_trial_accelerator; /**< \brief accelerator of the trial field */
};

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::result_t
	double_integral<Kernel, Test, Trial>::m_result;

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::test_accelerator_t
	double_integral<Kernel, Test, Trial>::m_test_accelerator;

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::trial_accelerator_t
	double_integral<Kernel, Test, Trial>::m_trial_accelerator;

#endif

