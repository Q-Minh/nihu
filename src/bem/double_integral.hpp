/**
* \file double_integral.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class double_integral and its specialisations
*/
#ifndef DOUBLE_INTEGRAL_HPP_INCLUDED
#define DOUBLE_INTEGRAL_HPP_INCLUDED

#include "../util/plain_type.hpp"
#include "kernel.hpp"
#include "field_type_accelerator.hpp"
#include "singular_accelerator.hpp"

/**
 * \brief class evaluating double integrals of the weighted residual approach
 * \tparam Kernel type of the kernel to integrate
 * \tparam Test type of the test field
 * \tparam Trial type of the trial field
 */
template <class Kernel, class Test, class Trial>
class double_integral
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"Kernel must be derived from kernel_base<Kernel>");
public:
	typedef Kernel kernel_t;		/**< \brief template parameter as nested type */
	typedef Test test_field_t;		/**< \brief template parameter as nested type */
	typedef Trial trial_field_t;	/**< \brief template parameter as nested type */

	typedef typename kernel_t::input_t kernel_input_t;		/**< \brief input type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;	/**< \brief result type of kernel */

	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;

	typedef field_type_accelerator<test_field_t, quadrature_family_t> test_field_type_accelerator_t;	/**< \brief type of the accelerator of the test field */
	typedef field_type_accelerator_pool<test_field_t, quadrature_family_t> test_field_type_accelerator_pool_t;	/**< \brief type of the accelerator pool of the test field */
	typedef field_type_accelerator<trial_field_t, quadrature_family_t> trial_field_type_accelerator_t;	/**< \brief type of the accelerator of the trial field */
	typedef field_type_accelerator_pool<trial_field_t, quadrature_family_t> trial_field_type_accelerator_pool_t;	/**< \brief type of the accelerator pool of the trial field */

	typedef singular_accelerator<kernel_t, test_field_t, trial_field_t> singular_accelerator_t;

	typedef typename test_field_t::nset_t::shape_t test_shape_t;	/**< \brief type of test shape function */
	typedef typename trial_field_t::nset_t::shape_t trial_shape_t;	/**< \brief type of trial shape function */

	typedef typename plain_type<
		typename product_type<
			kernel_result_t,
			typename product_type<test_shape_t, Eigen::Transpose<trial_shape_t> >::type
		>::type
	>::type result_t;

	/** \brief evaluate double integral with selected accelerators
	 * \param [in] test_field the test field to integrate on
	 * \param [in] test_acc field type accelerator of the test field
	 * \param [in] trial_field the trial field to integrate on
	 * \param [in] trial_acc field type accelerator of the trial field
	 * \return reference to the integration result
	 */
	static result_t const &eval_on_accelerator(
		test_field_t const &test_field,
		test_field_type_accelerator_t const &test_acc,
		trial_field_t const &trial_field,
		trial_field_type_accelerator_t const &trial_acc)
	{
		for(auto test_it = test_acc.cbegin(); test_it != test_acc.cend(); ++test_it)
		{
			kernel_input_t test_input(test_field.get_elem(), test_it->get_quadrature_elem());
			auto left = (test_it->get_shape() * test_input.get_jacobian()).eval();
			for(auto trial_it = trial_acc.cbegin(); trial_it != trial_acc.cend(); ++trial_it)
			{
				kernel_input_t trial_input(trial_field.get_elem(), trial_it->get_quadrature_elem());
				auto right = (trial_it->get_shape().transpose() * trial_input.get_jacobian()).eval();
				m_result += left * kernel_t::eval(test_input, trial_input) * right;
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
		m_result = result_t();	// clear result

		// check singularity
		if (m_singular_accelerator.eval(test_field, trial_field) != REGULAR)
		{
			// TODO: calculate
			return m_result;
		}

		// select quadrature
		kernel_input_t test_center(test_field.get_elem(), m_test_field_accelerator_pool[0]->cbegin()->get_quadrature_elem());
		kernel_input_t trial_center(trial_field.get_elem(), m_trial_field_accelerator_pool[0]->cbegin()->get_quadrature_elem());
		unsigned degree = kernel_t::estimate_complexity(test_center, trial_center);

		return eval_on_quadrature(test_field, *(m_test_field_accelerator_pool[degree]), trial_field, *(m_trial_field_accelerator_pool[degree]));
	}

protected:
	static result_t m_result; /**< \brief the integral result stored as static variable */

	static const test_field_type_accelerator_pool_t m_test_field_accelerator_pool;
	static const trial_field_type_accelerator_pool_t m_trial_field_accelerator_pool;
	static const singular_accelerator_t m_singular_accelerator;
};

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::result_t
	double_integral<Kernel, Test, Trial>::m_result;

template<class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::test_field_type_accelerator_pool_t
	const double_integral<Kernel, Test, Trial>::m_test_field_accelerator_pool;

template<class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::trial_field_type_accelerator_pool_t
	const double_integral<Kernel, Test, Trial>::m_trial_field_accelerator_pool;

template<class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::singular_accelerator_t
	const double_integral<Kernel, Test, Trial>::m_singular_accelerator;


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
	typedef typename kernel_t::result_t kernel_result_t;		/**< \brief input type of kernel */

	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;		/**< \brief quadrature family type */

	typedef typename quadrature_domain_traits<quadrature_family_t, typename trial_field_t::elem_t::domain_t>::quadrature_type trial_quadrature_t;	/**< \brief type of trial quadrature */
	typedef typename trial_quadrature_t::quadrature_elem_t quadrature_elem_t;	/**< \brief type of quadrature element */

	typedef field_type_accelerator<trial_field_t, quadrature_family_t> trial_field_type_accelerator_t;	/**< \brief type of the accelerator of the trial field */
	typedef field_type_accelerator_pool<trial_field_t, quadrature_family_t> trial_field_type_accelerator_pool_t;	/**< \brief type of the accelerator pool of the trial field */

	typedef singular_accelerator<kernel_t, test_field_t, trial_field_t> singular_accelerator_t;

	typedef typename test_field_t::nset_t test_nset_t;		/**< \brief type of element's N-set */
	typedef typename trial_field_t::nset_t::shape_t trial_shape_t;		/**< \brief type of element's N-set */

	typedef typename plain_type<
		typename product_type<kernel_result_t, Eigen::Transpose<trial_shape_t> >::type
	>::type result_t;

	/** \brief evaluate double integral with selected trial field quadrature
	 * \param [in] test_field the test field to integrate on
	 * \param [in] trial_field the trial field to integrate on
	 * \param [in] trial_quad the trial accelerator
	 * \return reference to the integration result
	 */
	static result_t const &eval_on_accelerator(
		test_field_t const &test_field,
		trial_field_t const &trial_field,
		trial_field_type_accelerator_t const &trial_acc)
	{
		for(auto test_it = test_nset_t::corner_begin(); test_it != test_nset_t::corner_end(); ++test_it)
		{
			quadrature_elem_t qe(*test_it);
			kernel_input_t collocational_point(test_field.get_elem(), qe);
			for(auto trial_it = trial_acc.cbegin(); trial_it != trial_acc.cend(); ++trial_it)
			{
				kernel_input_t trial_input(trial_field.get_elem(), trial_it->get_quadrature_elem());

				m_result.row(test_it - test_nset_t::corner_begin()) +=
					kernel_t::eval(collocational_point, trial_input) *
					(trial_input.get_jacobian() * trial_it->get_shape().transpose());
			}
		}

		return m_result;
	}

	template <class Quadrature>
	static result_t const &eval_on_quadrature(
		test_field_t const &test_field,
		trial_field_t const &trial_field,
		Quadrature const &trial_quad)
	{
		for(auto test_it = test_nset_t::corner_begin(); test_it != test_nset_t::corner_end(); ++test_it)
		{
			quadrature_elem_t qe(*test_it);
			kernel_input_t collocational_point(test_field.get_elem(), qe);
			// Go through trial quadrature
			for(auto trial_it = trial_quad.begin(); trial_it != trial_quad.end(); ++trial_it)
			{
				// Create kernel input
				kernel_input_t trial_input(trial_field.get_elem(), *trial_it);

				m_result.row(test_it - test_nset_t::corner_begin()) +=
					kernel_t::eval(collocational_point, trial_input) *
					(trial_input.get_jacobian() *
					 trial_field_t::nset_t::eval_shape(trial_it->get_xi()).transpose());
					//trial_it->get_shape().transpose());
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
		m_result = result_t();	// clear result

		// check singularity
		if (m_singular_accelerator.eval(test_field, trial_field) != REGULAR)
			return eval_on_quadrature(
				test_field,
				trial_field,
				*m_singular_accelerator.get_trial_quadrature()
			);

		// select quadrature
		quadrature_elem_t qe(test_field_t::elem_t::domain_t::get_center());
		kernel_input_t test_center(test_field.get_elem(), qe);
		kernel_input_t trial_center(trial_field.get_elem(), m_trial_field_accelerator_pool[0]->cbegin()->get_quadrature_elem());
		unsigned degree = kernel_t::estimate_complexity(test_center, trial_center);
		return eval_on_accelerator(test_field, trial_field, *(m_trial_field_accelerator_pool[degree]));
	}

protected:
	static result_t m_result; /**< \brief the integral result stored as static variable */
	static const trial_field_type_accelerator_pool_t m_trial_field_accelerator_pool;
	static const singular_accelerator_t m_singular_accelerator;
};

template<class Kernel, class Trial, class ElemType, class FieldOption>
typename double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::result_t
	double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::m_result;

template<class Kernel, class Trial, class ElemType, class FieldOption>
typename double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::trial_field_type_accelerator_pool_t
	const double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::m_trial_field_accelerator_pool;

template<class Kernel, class Trial, class ElemType, class FieldOption>
typename double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::singular_accelerator_t
	const double_integral<Kernel, field<ElemType, FieldOption, dirac_field>, Trial>::m_singular_accelerator;

#endif

