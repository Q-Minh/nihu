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
#include "elem_accelerator.hpp"
#include "function_space.hpp"

/**
 * \brief integrates a kernel over a field and stores the result in a static variable
 * \details stores the result in a static variable
 * \tparam Field the field type that needs to be handled
 * \tparam Kernel the kernel type that needs to be evaluated
 */
template <class Field, class Kernel>
class weighted_field_integral
{
public:
	/** \brief template parameter as nested type */
	typedef Field field_t;
	/** \brief template parameter as nested type */
	typedef Kernel kernel_t;

	/** \brief the input type of the kernel */
	typedef typename kernel_t::input_t kernel_input_t;
	/** \brief the result type of the kernel */
	typedef typename kernel_t::result_t kernel_result_t;

	/** \brief type of the quadrature */
	typedef gauss_quadrature<typename field_t::elem_t::domain_t> quadrature_t;
	/** \brief type of a quadrature element */
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;

	/** \brief type of the element's N-set */
	typedef typename field_t::nset_t nset_t;
	/** \brief type of the weighted kernel result */
	typedef typename Eigen::Matrix<typename kernel_t::scalar_t, field_t::num_dofs, kernel_t::num_elements> result_t;

	/**
	 * \brief evaluate the integral over a specific field with a quadrature
	 * \param field the field over which integration is performed
	 * \param q the quadrature
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	static result_t const &eval(field_t const &field, quadrature_t const &q)
	{
		// clear result
		result = result_t();
		std::for_each(q.begin(), q.end(),
			[&] (quadrature_elem_t const &qe) {
			// construct kernel input
			kernel_input_t input(field.get_elem(), qe);
			// evaluate kernel and get reference to result
			auto kernel_res = kernel_t::eval(input);
			// multiply kernel result with weighting matrix and quadrature weight
			result += nset_t::eval_L(qe.get_xi()) * (kernel_res * (input.get_jacobian()));
		});
		return result;
	}

protected:
	/** \brief the integral result stored as static variable */
	static result_t result;
};

/** \brief definition of the static integral result */
template <class Field, class Kernel>
typename weighted_field_integral<Field, Kernel>::result_t
	weighted_field_integral<Field, Kernel>::result;


/**
 * \brief  integrates a kernel over a function space and stores the result in static variable
 * \tparam FunctionSpace the function space over which integration is performed
 * \tparam Kernel the kernel to be integrated
 */
template <class FunctionSpace, class Kernel>
class weighted_integral
{
public:
	/** \brief template paramter as nested type */
	typedef FunctionSpace function_space_t;
	/** \brief template paramter as nested type */
	typedef Kernel kernel_t;

	/** \brief the class type abbreviated */
	typedef weighted_integral<function_space_t, kernel_t> weighted_integral_t;

	/** \brief the kernel input type */
	typedef typename kernel_t::input_t kernel_input_t;
	typedef elem_accelerator<kernel_input_t> elem_accelerator_t;

	/** \brief the element type vector inherited from the function space */
	typedef typename function_space_t::elem_type_vector_t elem_type_vector_t;
	/** \brief integration result type */
	typedef Eigen::Matrix<typename kernel_t::scalar_t, Eigen::Dynamic, kernel_t::num_elements> result_vector_t;

protected:
	/**
	 * \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class elem_t>
	struct eval_on_elemtype
	{
		// internal struct type is needed because tmp::call_each works on metafunctions
		struct type
		{
			/** \brief the field type */
			typedef field<elem_t, typename function_space_t::field_option> field_t;
			/** \brief the field integrator type */
			typedef weighted_field_integral<field_t, kernel_t> weighted_field_integral_t;
			/** \brief result type of field integrator */
			typedef typename weighted_field_integral_t::result_t result_t;
			/** \brief dofs type needed to iterate */
			typedef typename field_t::dofs_t dofs_t;

			/** \brief type of the quadrature */
			typedef gauss_quadrature<typename elem_t::domain_t> quadrature_t;

			void operator() (weighted_integral_t &wi)
			{
				// create quadrature (expensive, memory is allocated and eigenvalue problem is solved)
				quadrature_t q(2);
				
				// integrate for each element of the same type
				std::for_each(
					wi.func_space.template begin<elem_t>(),
					wi.func_space.template end<elem_t>(),
					[&wi, &q] (field_t const &f)
				{
					// get reference to field integral result
					result_t const &I = weighted_field_integral_t::eval(f, q);
					// write result into result vector
					dofs_t const &dofs = f.get_dofs();
					for (unsigned i = 0; i < field_t::num_dofs; ++i)
						wi.result_vector.row(dofs(i)) += I.row(i);
				}
				);
			}
		};
	};

public:
	/**
	 * \brief constructor initialises the function space reference member and allocates for the result
	 * \param func_space the function space over which integration is performed
	 */
	weighted_integral(function_space_t const &func_space) : func_space(func_space)
	{
		result_vector.resize(func_space.get_num_dofs(), Eigen::NoChange);
	}

	/** \brief this works only for quad_1_elem elements */
	void accelerate(void)
	{
		typedef quad_1_elem elem_t;
		typedef gauss_quadrature<elem_t::domain_t> quadrature_t;
		typedef field<elem_t,  typename function_space_t::field_option> field_t;

		// create quadrature (expensive, memory is allocated and eigenvalue problem is solved)
		quadrature_t q(2);
				
		// accelerate for each element of the same type
		std::for_each(
			func_space.template begin<elem_t>(),
			func_space.template end<elem_t>(),
			[this, &q] (field_t const &f)
		{
				elem_acc.push_back(
					elem_accelerator_t(f.get_elem(), q)
					);
		}
		);
	}

	/**
	 * \brief evaluate integral and return reference to result vector
	 * \return reference to static member result
	 */
	result_vector_t const &eval(void)
	{
		result_vector = result_vector_t::Zero(result_vector.rows(), result_vector.ColsAtCompileTime);
		tmp::call_each<
			elem_type_vector_t,
			eval_on_elemtype<tmp::_1>,
			weighted_integral_t &
		>(*this);
		return result_vector;
	}

protected:
	/** \brief reference to the function space */
	function_space_t const &func_space;
	/** \brief the result vector */
	result_vector_t result_vector;

	EIGENSTDVECTOR(elem_accelerator_t) elem_acc;
};

#endif

