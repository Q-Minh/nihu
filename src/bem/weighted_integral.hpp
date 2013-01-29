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
	typedef typename Eigen::Matrix<				/**< \brief type of weighted kernel result */
		typename kernel_t::scalar_t, field_t::num_dofs, kernel_t::num_elements
	> result_t;

	/**
	 * \brief evaluate the integral over a specific field with a quadrature
	 * \param field the field over which integration is performed
	 * \param q the quadrature
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	static result_t const &eval(field_t const &field, quadrature_t const &q)
	{
		result = result_t();
		for(auto it = q.begin(); it != q.end(); ++it)
		{
			kernel_input_t input(field.get_elem(), *it);
			auto kernel_res = kernel_t::eval(input);
			result += nset_t::eval_L(it->get_xi()) * (kernel_res * (input.get_jacobian()));
		}
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
 * \brief integrates a kernel over a field and stores the result in a static variable
 * \tparam Field the field type that needs to be handled
 * \tparam Kernel the kernel type that needs to be evaluated
 */
template <class Field, class Kernel>
class weighted_accelerated_field_integral : public weighted_field_integral<Field, Kernel>
{
public:
	typedef weighted_field_integral<Field, Kernel> base;

	typedef typename base::result_t result_t;
	typedef typename base::kernel_t kernel_t;
	typedef typename base::kernel_input_t kernel_input_t;
	using base::result;

	/**
	 * \brief evaluate the integral over a specific field with a quadrature
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	template <class kernel_iter_t, class nset_iter_t>
	static result_t const &eval(kernel_iter_t k_begin, kernel_iter_t k_end, nset_iter_t n_begin)
	{
		for (result = result_t(); k_begin != k_end; ++k_begin, ++n_begin)
			result += (*n_begin) * (kernel_t::eval(*k_begin) * k_begin->get_jacobian());
		return result;
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
	/** \brief template paramter as nested type */
	typedef FunctionSpace function_space_t;
	/** \brief template paramter as nested type */
	typedef Kernel kernel_t;

	/** \brief the class type abbreviated */
	typedef weighted_integral<function_space_t, kernel_t> weighted_integral_t;

	typedef typename function_space_t::field_option_t field_option_t;

	/** \brief the kernel input type */
	typedef typename kernel_t::input_t kernel_input_t;
	typedef elem_accelerator<kernel_input_t> elem_accelerator_t;

	/** \brief the element type vector inherited from the function space */
	typedef typename function_space_t::elem_type_vector_t elem_type_vector_t;
	/** \brief integration result type */
	typedef Eigen::Matrix<typename kernel_t::scalar_t, Eigen::Dynamic, kernel_t::num_elements> result_vector_t;

	typedef accelerator<elem_type_vector_t, field_option_t, kernel_input_t> accelerator_t;

protected:
	/**
	 * \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class elem_t>
	struct eval_on
	{
		// internal struct type is needed because tmp::call_each works on metafunctions
		struct type
		{
			/** \brief the field type */
			typedef field<elem_t, typename function_space_t::field_option_t> field_t;
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

	/**
	 * \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class elem_t>
	struct accelerated_eval_on
	{
		// internal struct type is needed because tmp::call_each works on metafunctions
		struct type : public eval_on<elem_t>::type
		{
			typedef typename eval_on<elem_t>::type base;

			typedef typename base::field_t field_t;
			typedef typename base::dofs_t dofs_t;
			typedef typename base::result_t result_t;

			/** \brief the field integrator type */
			typedef weighted_accelerated_field_integral<field_t, kernel_t> weighted_field_integral_t;

			void operator() (weighted_integral_t &wi)
			{
				auto nset_it = wi.accel.template get_nset_pool<elem_t>()[1].begin();

				typename function_space_t::template field_iterator_t<elem_t> field_it = wi.func_space.template begin<elem_t>();
				for (auto acc_it = wi.accel.template elem_begin<elem_t>();
					acc_it != wi.accel.template elem_end<elem_t>();
					++acc_it, ++field_it)
				{
					result_t const &I = weighted_field_integral_t::eval((*acc_it)[1].begin(), (*acc_it)[1].end(), nset_it);
					// write result into result vector
					dofs_t const &dofs = (*field_it).get_dofs();
					for (unsigned i = 0; i < field_t::num_dofs; ++i)
						wi.result_vector.row(dofs(i)) += I.row(i);
				}
			}
		};
	};


	template <class elem_t>
	struct accelerate_on
	{
		struct type
		{
			void operator() (weighted_integral_t &wi)
			{
				std::for_each(
					wi.func_space.template begin<elem_t>(),
					wi.func_space.template end<elem_t>(),
					[&] (field<elem_t, typename function_space_t::field_option_t> const &f) {
						wi.accel.add_elem(f.get_elem());
				});
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
		accelerate();
	}

	void accelerate(void)
	{
		tmp::call_each<
			elem_type_vector_t,
			accelerate_on<tmp::_1>,
			weighted_integral_t &
		>(*this);
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
			accelerated_eval_on<tmp::_1>,
			weighted_integral_t &
		>(*this);
		return result_vector;
	}

protected:
	
	function_space_t const &func_space;	/**< \brief reference to the function space */
	result_vector_t result_vector;		/**< \brief the result vector */

	accelerator_t accel;
};

#endif


