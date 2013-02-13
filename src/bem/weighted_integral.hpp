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

	/**
	 * \brief evaluate the regular integral over a specific field with a quadrature
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
			m_result += nset_t::eval_L(it->get_xi()) * (kernel_t::eval(input) * (input.get_jacobian()));
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


/**
 * \brief integrates a kernel over a pair of fields and stores the result in a static variable
 * \tparam Field1 the first field type that needs to be handled
 * \tparam Field2 the first field type that needs to be handled
 * \tparam Kernel the kernel type that needs to be evaluated
 */
template <class Field1, class Field2, class Kernel>
class weighted_double_field_integral
{
public:
	typedef Field1 field1_t;	/**< \brief template parameter as nested type */
	typedef Field2 field2_t;	/**< \brief template parameter as nested type */
	typedef Kernel kernel_t;/**< \brief template parameter as nested type */

	typedef typename kernel_t::input_t kernel_input_t;	/**< \brief input type of kernel */
	typedef typename kernel_t::result_t kernel_result_t;/**< \brief result type of kernel */

	typedef gauss_quadrature<typename field1_t::elem_t::domain_t> quadrature1_t; /**< \brief type of quadrature */
	typedef gauss_quadrature<typename field2_t::elem_t::domain_t> quadrature2_t; /**< \brief type of quadrature */
	typedef typename quadrature1_t::quadrature_elem_t quadrature1_elem_t;	/**< \brief type of quadrature element */
	typedef typename quadrature2_t::quadrature_elem_t quadrature2_elem_t;	/**< \brief type of quadrature element */

	typedef typename field1_t::nset_t nset1_t;	/**< \brief type of element's N-set */
	typedef typename field2_t::nset_t nset2_t;	/**< \brief type of element's N-set */
	typedef typename Eigen::Matrix<
		typename kernel_t::scalar_t, field1_t::num_dofs * field2_t::num_dofs, kernel_t::num_elements
	> result_t;									/**< \brief integration result type */
	typedef typename Eigen::Matrix<
		typename kernel_t::scalar_t, field2_t::num_dofs, kernel_t::num_elements
	> local_result_t;									/**< \brief local integration result type */

	/**
	 * \brief evaluate the integral over a specific field pair with a quadrature
	 * \param field1 the field over which integration is performed
	 * \param q_pool1 the quadrature pool
	 * \param field2 the field over which integration is performed
	 * \param q_pool2 the quadrature pool
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	template <class quad_pool1, class quad_pool2>
	static result_t const &eval(field1_t const &field1, quad_pool1 const &q_pool1,
								field2_t const &field2, quad_pool2 const &q_pool2)
	{
		m_result = result_t();	// clear result
		// select quadrature
		kernel_input_t trial1(field1.get_elem(), q_pool1[0][0]);
		kernel_input_t trial2(field2.get_elem(), q_pool2[0][0]);
		unsigned degree = kernel_t::estimate_complexity(trial1, trial2);

		for(auto it1 = q_pool1[degree].begin(); it1 != q_pool1[degree].end(); ++it1)
		{
			local_result_t result();

			kernel_input_t input1(field1.get_elem(), *it1);
			for(auto it2 = q_pool2[degree].begin(); it2 != q_pool2[degree].end(); ++it2)
			{
				kernel_input_t input2(field2.get_elem(), *it2);
				result += nset2_t::eval_L(it2->get_xi()) * (kernel_t::eval(input1, input2) * (input2.get_jacobian()));
			}
			result *= input1.get_jacobian();

			auto N1 = nset1_t::eval_L(it1->get_xi());
			for (int i = 0; i < field1_t::num_dofs; ++i)
				m_result.block<field2_t::num_dofs, kernel_t::num_elements>(i*field2_t::num_dofs,0) += result * N1[i];
		}
		return m_result;
	}

protected:
	static result_t m_result; /**< \brief the integral result stored colun-wise as static variable */
};

/** \brief definition of the static integral result */
template <class Field1, class Field2, class Kernel>
typename weighted_double_field_integral<Field1, Field2, Kernel>::result_t
	weighted_double_field_integral<Field1, Field2, Kernel>::m_result;




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
	struct integrate_on { struct type {
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
				result_t const &I = weighted_field_integral_t::eval(f, quadrature_pool);
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
	struct accelerated_integrate_on { struct type : public integrate_on<elem_t>::type
	{
		typedef typename integrate_on<elem_t>::type base;	/**< \brief base class abbreviation */

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


