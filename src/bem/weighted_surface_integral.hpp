/**
* \file weighted_surface_integral.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class weighted_surface_integral
*/
#ifndef WEIGHTED_SURFACE_INTEGRAL_HPP_INCLUDED
#define WEIGHTED_SURFACE_INTEGRAL_HPP_INCLUDED

#include <Eigen/Dense>

#include <vector>

#include "kernel.hpp"
#include "quadrature.hpp"
#include "descriptor.hpp"
#include "function_space.hpp"
#include "accelerator.hpp"

/** \brief integrates a kernel over a pair of fields
* \tparam Field1 the type of the first field
* \tparam Field2 the type of the second field
* \tparam Kernel the kernel type that needs to be evaluated
*/
template <class Kernel, class Field1, class Field2>
class weighted_field_surface_integral
{
public:
	typedef Field1 field1_t;	/**< \brief template parameter as nested type */
	typedef Field2 field2_t;	/**< \brief template parameter as nested type */
	typedef Kernel kernel_t;	/**< \brief template parameter as nested type */

	typedef typename kernel_t::input_t kernel_input_t;	/**< \brief input type of kernel */

	typedef gauss_quadrature<typename field1_t::elem_t::domain_t> quadrature1_t; /**< \brief type of quadrature */
	typedef gauss_quadrature<typename field2_t::elem_t::domain_t> quadrature2_t; /**< \brief type of quadrature */

	typedef typename field1_t::nset_t nset1_t;	/**< \brief type of element's N-set */
	typedef typename field2_t::nset_t nset2_t;	/**< \brief type of element's N-set */
	static const unsigned num_dofs1 = field1_t::num_dofs;
	static const unsigned num_dofs2 = field2_t::num_dofs;
	typedef typename field1_t::dofs_t dofs1_t;
	typedef typename field2_t::dofs_t dofs2_t;
	typedef typename Eigen::Matrix<
		typename kernel_t::scalar_t, num_dofs2, kernel_t::num_elements
	> inner_result_t;	/**< \brief integration result type */
	typedef typename nset1_t::shape_t N1_t;
	typedef typename nset2_t::shape_t N2_t;

	/** \brief evaluate the regular integral over a specific field pair
	* \tparam quad1_pool
	* \tparam quad2_pool
	* \tparam result_matrix_t
	* \param field1
	* \param q1_pool
	* \param field2
	* \param q2_pool
	* \param result_matrix 
	*/
	template <class quad1_pool, class quad2_pool, class result_matrix_t>
	static void eval(field1_t const &field1, quad1_pool const &q1_pool,
		field2_t const &field2, quad2_pool const &q2_pool,
		result_matrix_t &result_matrix)
	{
		dofs1_t const &dofs1 = field1.get_dofs();
		dofs2_t const &dofs2 = field2.get_dofs();

		// select quadrature
		kernel_input_t trial1(field1.get_elem(), q1_pool[0][0]);
		kernel_input_t trial2(field2.get_elem(), q2_pool[0][0]);
		unsigned degree = kernel_t::estimate_complexity(trial1, trial2);

		for(auto it1 = q1_pool[degree].begin(); it1 != q1_pool[degree].end(); ++it1)
		{
			inner_result_t inner_result; // clear inner result

			kernel_input_t input1(field1.get_elem(), *it1);
			for(auto it2 = q2_pool[degree].begin(); it2 != q2_pool[degree].end(); ++it2)
			{
				kernel_input_t input2(field2.get_elem(), *it2);
				N2_t const &N2 = nset2_t::eval_shape(it2->get_xi());
				double jac2 = input2.get_jacobian();
				inner_result += N2 * (kernel_t::eval(input1, input2) * jac2);
			}

			inner_result *= input1.get_jacobian();

			N1_t const &N1 = nset1_t::eval_shape(it1->get_xi());
			for (unsigned i = 0; i < num_dofs1; ++i)
				for (unsigned j = 0; j <num_dofs2; ++j)
					result_matrix(dofs1(i),dofs2(j)) += N1(i) * inner_result(j);
		}
	}
};

template <class field1_t, class field2_t>
class singularity_check
{
public:
	static bool eval(field1_t const &f1, field2_t const &f2)
	{
		return false;
	}
};

template <class field1_t>
class singularity_check<field1_t, field1_t>
{
public:
	static bool eval(field1_t const &f1, field1_t const &f2)
	{
		return &(f1.get_elem()) == &(f2.get_elem());
	}
};



/** \brief  integrates a kernel over a function space and stores the result in static variable
* \tparam FunctionSpace the function space over which integration is performed
* \tparam Kernel the kernel to be integrated
*/
template <class Kernel, class FunctionSpace1, class FunctionSpace2>
class weighted_surface_integral
{
public:
	typedef Kernel kernel_t;	/**< \brief template paramter as nested type */
	typedef FunctionSpace1 function_space1_t;	/**< \brief template paramter as nested type */
	typedef FunctionSpace2 function_space2_t;	/**< \brief template paramter as nested type */

	typedef typename function_space1_t::field_option_t field_option1_t;	/**< \brief field option parameter */
	typedef typename function_space2_t::field_option_t field_option2_t;	/**< \brief field option parameter */
	typedef typename kernel_t::input_t kernel_input_t;	/**< \brief kernel input type */
	typedef typename function_space1_t::elem_type_vector_t elem_type_vector_t;	/**< \brief the element type vector */

	typedef accelerator<elem_type_vector_t> accelerator_t;	/**< \brief type of accelerator class */

protected:
	template <class elem1_t, class elem2_t, class result_matrix_t>
	struct integrate_on { struct type {
		typedef field<elem1_t, field_option1_t> field1_t;		/**< \brief the field type */
		typedef field<elem2_t, field_option2_t> field2_t;		/**< \brief the field type */
		typedef weighted_surface_integral<kernel_t, function_space1_t, function_space2_t> weighted_surface_integral_t;	/**< \brief the class type abbreviated */
		typedef weighted_field_surface_integral<kernel_t, field1_t, field2_t> weighted_field_surface_integral_t;	/**< \brief the field integrator type */

		typedef gauss_quadrature<typename elem1_t::domain_t> quadrature1_t; /**< \brief type of quadrature */
		typedef gauss_quadrature<typename elem2_t::domain_t> quadrature2_t; /**< \brief type of quadrature */

		typedef singularity_check<field1_t, field2_t> is_singular;

		void operator() (weighted_surface_integral_t &wi, result_matrix_t &result_matrix)
		{
			auto quadrature1_pool = wi.m_accel.template get_quadrature_pool<elem1_t>();
			auto quadrature2_pool = wi.m_accel.template get_quadrature_pool<elem2_t>();

			for (auto f1 = wi.m_func_space1.template elem_begin<elem1_t>();
				f1 != wi.m_func_space1.template elem_end<elem1_t>();
				++f1)
			{

				for (auto f2 = wi.m_func_space2.template elem_begin<elem2_t>();
					f2 != wi.m_func_space2.template elem_end<elem2_t>();
					++f2)
				{
					if (is_singular::eval(*f1, *f2))
						continue;
					weighted_field_surface_integral_t::eval(*f1, quadrature1_pool, *f2, quadrature2_pool, result_matrix);
				}
			}
		}
	};};

public:
	weighted_surface_integral(function_space1_t const &func_space1,
		function_space2_t const &func_space2)
		: m_func_space1(func_space1),  m_func_space2(func_space2)
	{
	}

	/** \brief evaluate integral and return reference to result vector
	* \return reference to static member result
	*/
	template <class result_matrix_t>
	void integrate(result_matrix_t &result_matrix)
	{
		typename integrate_on<quad_1_elem, quad_1_elem, result_matrix_t>::type t;
		t(*this, result_matrix);
	}

protected:
	function_space1_t const &m_func_space1;	/**< \brief reference to the function space */
	function_space2_t const &m_func_space2;	/**< \brief reference to the function space */

	accelerator_t m_accel;	/**< \brief accelerator structure */
};

#endif
