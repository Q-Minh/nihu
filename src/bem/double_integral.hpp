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


template <bool isCollocational, class Kernel, class TestField, class TrialField>
struct accel_store
{
	typedef singular_accelerator<isCollocational, Kernel, TestField, TrialField> singular_accelerator_t;
	static singular_accelerator_t m_singular_accelerator;
};

template<bool isCollocational, class Kernel, class Test, class Trial>
typename accel_store<isCollocational, Kernel, Test, Trial>::singular_accelerator_t
	accel_store<isCollocational, Kernel, Test, Trial>::m_singular_accelerator;


template <class Field, class Family>
struct regular_pool_store
{
	typedef field_type_accelerator_pool<Field, Family> pool_t;
	static pool_t m_regular_pool;
};

template <class Field, class Family>
typename regular_pool_store<Field, Family>::pool_t
	regular_pool_store<Field, Family>::m_regular_pool;

/**
* \brief class evaluating double integrals of the weighted residual approach
* \tparam Kernel type of the kernel to integrate
* \tparam Test type of the test field
* \tparam Trial type of the trial field
*/
template <class Kernel, class TestField, class TrialField>
class double_integral
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"Kernel must be derived from kernel_base<Kernel>");
	static_assert(std::is_base_of<field_base<TestField>, TestField>::value,
		"TestField must be derived from field_base<TestField>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"TrialField must be derived from field_base<TrialField>");
public:
	typedef Kernel kernel_t;		/**< \brief template parameter as nested type */
	typedef TestField test_field_t;		/**< \brief template parameter as nested type */
	typedef TrialField trial_field_t;	/**< \brief template parameter as nested type */

	typedef typename kernel_t::input_t kernel_input_t;	/**< \brief input type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;	/**< \brief result type of kernel */

	/** \brief the quadrature family the kernel requires */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
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

	/** \brief result type of the weighted residual */
	typedef typename plain_type<
		typename product_type<
		kernel_result_t,
		typename product_type<test_shape_t, Eigen::Transpose<trial_shape_t> >::type
		>::type
	>::type result_t;

protected:
	/** \brief evaluate regular double integral with selected accelerators
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
		for (auto test_it = test_acc.cbegin(); test_it != test_acc.cend(); ++test_it)
		{
			kernel_input_t test_input(test_field.get_elem(), test_it->get_quadrature_elem());
			auto left = (test_it->get_shape() * test_input.get_jacobian()).eval();
			for (auto trial_it = trial_acc.cbegin(); trial_it != trial_acc.cend(); ++trial_it)
			{
				kernel_input_t trial_input(trial_field.get_elem(), trial_it->get_quadrature_elem());
				auto right = (trial_it->get_shape().transpose() * trial_input.get_jacobian()).eval();
				m_result += left * kernel_t::eval(test_input, trial_input) * right;
			}
		}

		return m_result;
	}

	/** \brief evaluate regular collocational integral with selected trial field accelerator
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \param [in] trial_acc the trial accelerator
	* \return reference to the integration result
	*/
	static result_t const &eval_collocational_on_accelerator(
		test_field_t const &test_field,
		trial_field_t const &trial_field,
		trial_field_type_accelerator_t const &trial_acc)
	{
		int row = 0;
		for (auto test_it = test_nset_t::corner_begin(); test_it != test_nset_t::corner_end(); ++test_it)
		{
			quadrature_elem_t qe(*test_it);	// quadrature weight is scalar_t(), but we don't care
			/** \todo we should not compute the Jacobian here */
			kernel_input_t collocational_point(test_field.get_elem(), qe);
			for (auto trial_it = trial_acc.cbegin(); trial_it != trial_acc.cend(); ++trial_it)
			{
				kernel_input_t trial_input(trial_field.get_elem(), trial_it->get_quadrature_elem());
				m_result.row(row) += kernel_t::eval(collocational_point, trial_input) *
					(trial_input.get_jacobian() * trial_it->get_shape().transpose());
			}
			++row;
		}

		return m_result;
	}


	/** \brief evaluate double singular integral with selected singular accelerator
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \param [in] begin begin iteartor of the singular quadrature
	* \param [in] end end iterator of the singular quadrature
	* \return reference to the integration result
	*/
	template <class singular_iterator_t>
	static result_t const &eval_singular_on_accelerator(
		test_field_t const &test_field,
		trial_field_t const &trial_field,
		singular_iterator_t begin,
		singular_iterator_t end)
	{
		while (begin != end)
		{
			kernel_input_t test_input(test_field.get_elem(), begin.get_test_quadrature_elem());
			kernel_input_t trial_input(trial_field.get_elem(), begin.get_trial_quadrature_elem());

			/** \todo check if lazy evaluation is still faster */
			auto left = (test_field_t::nset_t::eval_shape(begin.get_test_quadrature_elem().get_xi())
				* test_input.get_jacobian()).eval();
			auto right = (trial_field_t::nset_t::eval_shape(begin.get_trial_quadrature_elem().get_xi())
				* trial_input.get_jacobian()).eval();

			m_result += left * kernel_t::eval(test_input, trial_input) * right.transpose();

			++begin;
		}

		return m_result;
	}

	/** \brief evaluate collocational singular integral with selected singular accelerator
	* \tparam singular_iterator_t type of the singular quadrature iterator
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \param [in] begin the begin interator of the selected quadrature
	* \param [in] end end iterator of the selected quadrature
	* \return reference to the integration result
	*/
	template <class singular_accelerator_t>
	static result_t const &eval_collocational_singular_on_accelerator(
		test_field_t const &test_field,
		trial_field_t const &trial_field,
		singular_accelerator_t const &sa)
	{
		for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
		{
			quadrature_elem_t qe(test_nset_t::corner_at(idx));
			kernel_input_t collocational_point(test_field.get_elem(), qe);

			auto const &quad = sa.get_trial_quadrature(idx);

			for (auto quad_it = quad.begin(); quad_it != quad.end(); ++quad_it)
			{
				kernel_input_t trial_input(trial_field.get_elem(), *quad_it);

				m_result.row(idx) += kernel_t::eval(collocational_point, trial_input) *
					(trial_input.get_jacobian() *
					trial_field_t::nset_t::eval_shape(quad_it->get_xi()).transpose());
			}
		}

		return m_result;
	}

public:
	/** \brief evaluate double integral on given fields
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \return reference to the integration result
	*/
	static result_t const &eval(
		std::false_type,
		test_field_t const &test_field,
		trial_field_t const &trial_field)
	{
		typedef accel_store<false, kernel_t, test_field_t, trial_field_t> acc_store_t;
		auto &sa = acc_store_t::m_singular_accelerator;

		auto &test_ra = regular_test_store_t::m_regular_pool;
		auto &trial_ra = regular_trial_store_t::m_regular_pool;

		m_result.setZero();	// clear result

		// check singularity
		if (sa.is_singular(test_field, trial_field))
			return eval_singular_on_accelerator(
			test_field, trial_field, sa.begin(), sa.end());

		// select quadrature
		kernel_input_t test_center(test_field.get_elem(),
			test_ra[0]->cbegin()->get_quadrature_elem());
		kernel_input_t trial_center(trial_field.get_elem(),
			trial_ra[0]->cbegin()->get_quadrature_elem());

		unsigned degree = kernel_t::estimate_complexity(test_center, trial_center);

		return eval_on_accelerator(
			test_field, *(test_ra[degree]),
			trial_field, *(trial_ra[degree]));
	}


	/** \brief evaluate collocational integral on given fields
	* \param [in] test_field the test field to integrate on
	* \param [in] trial_field the trial field to integrate on
	* \return reference to the integration result
	*/
	static result_t const &eval(
		std::true_type,
		test_field_t const &test_field,
		trial_field_t const &trial_field)
	{
		typedef accel_store<true, kernel_t, test_field_t, trial_field_t> acc_store_t;
		typename acc_store_t::singular_accelerator_t &sa = acc_store_t::m_singular_accelerator;

		typedef regular_pool_store<trial_field_t, quadrature_family_t> regular_trial_store_t;
		auto &trial_ra = regular_trial_store_t::m_regular_pool;

		m_result.setZero();	// clear result

		// check singularity
		if (sa.is_singular(test_field, trial_field))
			return eval_collocational_singular_on_accelerator(
			test_field, trial_field, sa);

		// select quadrature
		quadrature_elem_t qe(test_field_t::elem_t::domain_t::get_center());
		kernel_input_t test_center(test_field.get_elem(), qe);
		kernel_input_t trial_center(trial_field.get_elem(),
			trial_ra[0]->cbegin()->get_quadrature_elem());

		unsigned degree = kernel_t::estimate_complexity(test_center, trial_center);

		return eval_collocational_on_accelerator(
			test_field, trial_field, *(trial_ra[degree]));
	}

protected:
	/** \brief the integral result stored as static variable */
	static result_t m_result;
};

template <class Kernel, class Test, class Trial>
typename double_integral<Kernel, Test, Trial>::result_t
	double_integral<Kernel, Test, Trial>::m_result;

#endif

