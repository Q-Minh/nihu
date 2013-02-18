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

/** \brief integrates a kernel over a field and stores the result in a static variable
 * \tparam Field the field type that needs to be handled
 * \tparam Kernel the kernel type that needs to be evaluated
 */
template <class Kernel, class Field>
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

	/** \brief evaluate the regular integral over a specific field with a quadrature
	 * \param field the field over which integration is performed
	 * \param q_pool the quadrature pool
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	template <class quad_pool>
	static result_t const &eval(field_t const &field, quad_pool const &q_pool)
	{
		m_result = result_t();	// clear result
		// select quadrature
		kernel_input_t trial(field.get_elem(), q_pool[0][0]);
		unsigned degree = kernel_t::estimate_complexity(trial);

		m_result = result_t();	// clear result
		for(auto it = q_pool[degree].begin(); it != q_pool[degree].end(); ++it)
		{
			kernel_input_t input(field.get_elem(), *it);
			m_result += nset_t::eval_shape(it->get_xi()) * (kernel_t::eval(input) * (input.get_jacobian()));
		}
		return m_result;
	}

protected:
	static result_t m_result; /**< \brief the integral result stored as static variable */
};

/** \brief definition of the static integral result */
template <class Field, class Kernel>
typename weighted_field_integral<Field, Kernel>::result_t
	weighted_field_integral<Field, Kernel>::m_result;


/** \brief  integrates a kernel over a function space and stores the result in static variable
 * \tparam FunctionSpace the function space over which integration is performed
 * \tparam Kernel the kernel to be integrated
 */
template <class Kernel, class FunctionSpace>
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

	typedef accelerator<elem_type_vector_t> accelerator_t;	/**< \brief type of accelerator class */

protected:
	/** \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class elem_t>
	struct integrate_on { struct type {
		typedef field<elem_t, typename function_space_t::field_option_t> field_t;		/**< \brief the field type */
		typedef weighted_field_integral<kernel_t, field_t> weighted_field_integral_t;	/**< \brief the field integrator type */
		typedef typename weighted_field_integral_t::result_t result_t;				/**< \brief result type of field integrator */
		typedef typename field_t::dofs_t dofs_t;										/**< \brief dofs type needed to iterate */

		typedef gauss_quadrature<typename elem_t::domain_t> quadrature_t; /**< \brief type of quadrature */

		/** \brief evaluate integral on one element type
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
				result_t const &I = weighted_field_integral_t::eval(f, quadrature_pool);
				// write result into result vector
				dofs_t const &dofs = f.get_dofs();
				for (unsigned i = 0; i < field_t::num_dofs; ++i)
					wi.m_result_vector.row(dofs(i)) += I.row(i);
			}
			);
		}
	};};

public:
	/** \brief constructor initialises the function space reference member and allocates the result
	 * \param func_space the function space over which integration is performed
	 */
	weighted_integral(function_space_t const &func_space) : m_func_space(func_space)
	{
		m_result_vector.resize(m_func_space.get_num_dofs(), Eigen::NoChange);
	}

	/** \brief evaluate integral and return reference to result vector
	 * \return reference to static member result
	 */
	result_vector_t const &integrate(void)
	{
		m_result_vector = result_vector_t::Zero(m_result_vector.rows(), m_result_vector.ColsAtCompileTime);

		tmp::call_each<
			elem_type_vector_t,
			integrate_on<tmp::_1>,
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


