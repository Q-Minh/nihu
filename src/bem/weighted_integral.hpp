/**
* \file weighted_integral.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class weighted_integral
*/
#ifndef WEIGHTED_INTEGRAL_HPP_INCLUDED
#define WEIGHTED_INTEGRAL_HPP_INCLUDED

#include <Eigen/Dense>

#include <vector>

#include "kernel.hpp"
#include "quadrature.hpp"
#include "descriptor.hpp"
#include "function_space.hpp"
#include "accelerator.hpp"

/**
 * \brief integrates a kernel over a field and stores the result in a static variable
 * \tparam Field the field type that needs to be handled
 * \tparam Kernel the kernel type that needs to be evaluated
 */
template <class Field, class Kernel>
class weighted_field_integral
{
public:
	typedef Field field_t;	/**< \brief template parameter as nested type */
	typedef Kernel kernel_t;/**< \brief template parameter as nested type */

	typedef typename kernel_t::input_t kernel_input_t;	/**< \brief input type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;/**< \brief result type of kernel */

	typedef gauss_quadrature<typename field_t::elem_t::domain_t> quadrature_t; /**< \brief type of quadrature */
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;	/**< \brief type of quadrature element */

	typedef typename field_t::nset_t nset_t;	/**< \brief type of element's N-set */
	typedef typename Eigen::Matrix<
		typename kernel_t::scalar_t, field_t::num_dofs, kernel_t::num_elements
	> result_t;									/**< \brief integration result type */

protected:
	/**
	 * \brief evaluate the regular integral over a specific field with a quadrature
	 * \param field the field over which integration is performed
	 * \param q_pool the quadrature pool
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	template <class quad_pool>
	static result_t const &eval_regular_impl(field_t const &field, quad_pool const &q_pool)
	{
		// select quadrature
		kernel_input_t trial(field.get_elem(), q_pool[0][0]);
		unsigned degree = kernel_t::estimate_complexity(trial);

		m_result = result_t();	// clear result
		for(auto it = q_pool[degree].begin(); it != q_pool[degree].end(); ++it)
		{
			kernel_input_t input(field.get_elem(), *it);
			m_result += nset_t::eval_L(it->get_xi()) * (kernel_t::eval(input) * (input.get_jacobian()));
		}
		return m_result;
	}

public:
	/**
	 * \brief evaluate the regular integral over a specific field with a quadrature
	 * \param field the field over which integration is performed
	 * \param q_pool the quadrature pool
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	template <class quad_pool>
	static result_t const &eval_regular(field_t const &field, quad_pool const &q_pool)
	{
		m_result = result_t();	// clear result
		return eval_regular_impl(field, q_pool);
	}

	/**
	 * \brief evaluate the possibly singular integral over a specific field with a quadrature
	 * \param field the field over which integration is performed
	 * \param source_dof the source dof
	 * \param q_pool the quadrature pool
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	template <class quad_pool>
	static result_t const &eval_surface(field_t const &field, unsigned source_dof, quad_pool const &q_pool)
	{
		m_result = result_t();	// clear result

		// check whether singular quadrature is needed
		typename field_t::dofs_t dofs = field.get_dofs();
		for (unsigned i = 0; i < field_t::num_dofs; ++i)
			if (dofs[i] == source_dof)
				return m_result;

		return eval_regular_impl(field, q_pool);
	}

protected:
	static result_t m_result; /**< \brief the integral result stored as static variable */
};

/** \brief definition of the static integral result */
template <class Field, class Kernel>
typename weighted_field_integral<Field, Kernel>::result_t
	weighted_field_integral<Field, Kernel>::m_result;

/**
 * \brief integrates a kernel over a field and stores the result in a static variable
 * \tparam Field the field type that needs to be handled
 * \tparam Kernel the kernel type that needs to be evaluated
 */
template <class Field, class Kernel>
class weighted_accelerated_field_integral : public weighted_field_integral<Field, Kernel>
{
public:
	typedef weighted_field_integral<Field, Kernel> base;	/**< \brief type of the base class */
	typedef typename base::result_t result_t;				/**< \brief integration result type */
	typedef typename base::kernel_t kernel_t;				/**< \brief kernel type */
	typedef typename base::kernel_input_t kernel_input_t;	/**< \brief kernel input type */

	/**
	 * \brief evaluate the integral over a specific field with a quadrature
	 * \param k_begin begin iterator of kernel input
	 * \param k_end end iterator of kernel input
	 * \param n_begin begin iterator of shape set values
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	template <class kernel_iter_t, class nset_iter_t>
	static result_t const &eval(kernel_iter_t k_begin, kernel_iter_t k_end, nset_iter_t n_begin)
	{
		base::m_result = result_t();
		for (;k_begin != k_end; ++k_begin, ++n_begin)
			base::m_result += (*n_begin) * (kernel_t::eval(*k_begin) * k_begin->get_jacobian());
		return base::m_result;
	}
};


/**
 * \brief  integrates a kernel over a function space and stores the result in static variable
 * \tparam FunctionSpace the function space over which integration is performed
 * \tparam Kernel the kernel to be integrated
 */
template <class FunctionSpace, class Kernel>
class weighted_integral
{
public:
	typedef FunctionSpace function_space_t;	/**< \brief template paramter as nested type */
	typedef Kernel kernel_t;	/**< \brief template paramter as nested type */

	typedef weighted_integral<function_space_t, kernel_t> weighted_integral_t;	/**< \brief the class type abbreviated */
	typedef typename function_space_t::field_option_t field_option_t;	/**< \brief field option parameter */
	typedef typename kernel_t::input_t kernel_input_t;	/**< \brief kernel input type */
	typedef typename function_space_t::elem_type_vector_t elem_type_vector_t;	/**< \brief the element type vector */
	typedef Eigen::Matrix<typename kernel_t::scalar_t, Eigen::Dynamic, kernel_t::num_elements> result_vector_t;	/**< \brief integration result type */

	typedef accelerator<elem_type_vector_t, field_option_t, kernel_input_t> accelerator_t;	/**< \brief type of accelerator class */

protected:
	/**
	 * \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class elem_t>
	struct eval_on { struct type {
		typedef field<elem_t, typename function_space_t::field_option_t> field_t;		/**< \brief the field type */
		typedef weighted_field_integral<field_t, kernel_t> weighted_field_integral_t;	/**< \brief the field integrator type */
		typedef typename weighted_field_integral_t::result_t result_t;				/**< \brief result type of field integrator */
		typedef typename field_t::dofs_t dofs_t;										/**< \brief dofs type needed to iterate */

		typedef gauss_quadrature<typename elem_t::domain_t> quadrature_t; /**< \brief type of quadrature */

		/**
		 * \brief evaluate integral on one element type
		 * \param wi constant reference to weighted integral object
		 */
		void operator() (weighted_integral_t &wi)
		{
			auto quadrature_pool = wi.m_accel.template get_quadrature_pool<elem_t>();
			// integrate for each element of the same type
			std::for_each(
				wi.m_func_space.template elem_begin<elem_t>(),
				wi.m_func_space.template elem_end<elem_t>(),
				[&] (field_t const &f)
			{
				// get reference to field integral result
				result_t const &I = weighted_field_integral_t::eval_regular(f, quadrature_pool);
				// write result into result vector
				dofs_t const &dofs = f.get_dofs();
				for (unsigned i = 0; i < field_t::num_dofs; ++i)
					wi.m_result_vector.row(dofs(i)) += I.row(i);
			}
			);
		}
	};};

	/**
	 * \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class elem_t>
	struct accelerated_eval_on { struct type : public eval_on<elem_t>::type
	{
		typedef typename eval_on<elem_t>::type base;	/**< \brief base class abbreviation */

		typedef typename base::field_t field_t;	/**< \brief type of elem-field */
		typedef typename base::dofs_t dofs_t;		/**< \brief type of DOF vector */
		typedef typename base::result_t result_t;	/**< \brief type of integration result */

		/** \brief the field integrator type */
		typedef weighted_accelerated_field_integral<field_t, kernel_t> weighted_field_integral_t;

		/**
		 * \brief evaluate integral on one element type with acceleration
		 * \param wi constant reference to weighted integral object
		 */
		void operator() (weighted_integral_t &wi)
		{
			auto nset_it = wi.m_accel.template get_nset_pool<elem_t>()[1].begin();

			auto field_it = wi.m_func_space.template begin<elem_t>();
			for (auto acc_it = wi.m_accel.template elem_begin<elem_t>();
				acc_it != wi.m_accel.template elem_end<elem_t>();
				++acc_it, ++field_it)
			{
				result_t const &I = weighted_field_integral_t::eval((*acc_it)[1].begin(), (*acc_it)[1].end(), nset_it);
				// write result into result vector
				dofs_t const &dofs = (*field_it).get_dofs();
				for (unsigned i = 0; i < field_t::num_dofs; ++i)
					wi.m_result_vector.row(dofs(i)) += I.row(i);
			}
		}
	};};


	template <class elem_t>
	struct accelerate_on {	struct type
	{
		/**
		 * \brief perform acceleration for one element type
		 * \param wi constant reference to weighted integral object
		 */
		void operator() (weighted_integral_t &wi)
		{
			for (auto it = wi.m_func_space.template begin<elem_t>();
					it != wi.m_func_space.template end<elem_t>();
					++it)
				wi.m_accel.add_elem((*it).get_elem());
		}
	};};

public:
	/**
	 * \brief constructor initialises the function space reference member and allocates the result
	 * \param func_space the function space over which integration is performed
	 */
	weighted_integral(function_space_t const &func_space) : m_func_space(func_space)
	{
		m_result_vector.resize(m_func_space.get_num_dofs(), Eigen::NoChange);

		/*
		tmp::call_each<
			elem_type_vector_t,
			accelerate_on<tmp::_1>,
			weighted_integral_t &
		>(*this);
		*/
	}

	/**
	 * \brief evaluate integral and return reference to result vector
	 * \return reference to static member result
	 */
	result_vector_t const &eval(void)
	{
		m_result_vector = result_vector_t::Zero(m_result_vector.rows(), m_result_vector.ColsAtCompileTime);

		tmp::call_each<
			elem_type_vector_t,
			eval_on<tmp::_1>,
			weighted_integral_t &
		>(*this);

		return m_result_vector;
	}

protected:
	function_space_t const &m_func_space;	/**< \brief reference to the function space */
	result_vector_t m_result_vector;	/**< \brief the result vector */

	accelerator_t m_accel;	/**< \brief accelerator structure */
};

#endif


