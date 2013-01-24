/**
* \file weighted_integral.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class weighted_integral
*/
#ifndef WEIGHTED_INTEGRAL_HPP_INCLUDED
#define WEIGHTED_INTEGRAL_HPP_INCLUDED

#include "kernel.hpp"
#include "quadrature.hpp"
#include "descriptor.hpp"
#include "function_space.hpp"

/**
 * \brief integrates a kernel over a field with a quadrature of specified order
 * \details stores the result in a static variable
 * \tparam Field the field type that needs to be handled
 * \tparam Kernel the kernel type that needs to be evaluated
 * \tparam N the order of quadrature
 */
template <class Field, class Kernel, unsigned N>
class weighted_field_integral
{
public:
	/** \brief template parameter as nested type */
	typedef Field field_t;
	/** \brief template parameter as nested type */
	typedef Kernel kernel_t;
	/** \brief template parameter as nested constant */
	static unsigned const quad_n = N;

	/** \brief the input type of the kernel */
	typedef typename kernel_t::input_t kernel_input_t;
	/** \brief the result type of the kernel */
	typedef typename kernel_t::result_t kernel_result_t;

	/** \brief type of the quadrature */
	typedef gauss_quad<typename field_t::elem_t::domain_t, N> quadrature_t;
	/** \brief type of the quadrature base points vector */
	typedef typename quadrature_t::xivec_t xivec_t;
	/** \brief type of the quadrature weights vector */
	typedef typename quadrature_t::weightvec_t weightvec_t;
	/** \brief type of the element's N-set */
	typedef typename field_t::nset_t nset_t;
	/** \brief type of the weighted kernel result */
	typedef typename Eigen::Matrix<typename kernel_t::scalar_t, field_t::num_dofs, 1> result_t;

	/**
	 * \brief evaluate the integral over a specific field
	 * \param field the field over which the integration is performed
	 * \return reference to the result of the integral. The actual result is stored in static variable
	 */
	static result_t const &eval(field_t const &field)
	{
		// get quadrature base points and weights in Eigen Matrix
		xivec_t xivec = quadrature_t::get_xi();
		weightvec_t weightvec = quadrature_t::get_weight();

		// clear result vector
		result = result_t();
		// for each quadrature point
		for (typename xivec_t::Index i = 0; i < quadrature_t::size; ++i)
		{
			// auto is used to have optimal performance (reference to a row)
			auto xi = xivec.row(i);
			// evaluate the kernel input and store the result
			kernel_input_t input(field.get_elem(), xi);
			// evaluate the kernel and store a reference to the result
			auto kernel_res = kernel_t::eval(input);
			// accumulate result multiplied by shape function weights and jacobian
			result += nset_t::eval_L(xi) * (kernel_res * (input.get_jacobian() * weightvec(i)));
		}

		// return reference to the stored result
		return result;
	}

protected:
	/** \brief the integral result stored as static variable */
	static result_t result;
};

/** \brief definition of the integral result */
template <class Field, class Kernel, unsigned N>
typename weighted_field_integral<Field, Kernel, N>::result_t
	weighted_field_integral<Field, Kernel, N>::result;


/**
 * \brief  integrates a kernel over a function space and stores the result in static variable
 * \tparam FunctionSpace the function space over which the integration is performed
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

	/** \brief the element type vector inherited from the function space */
	typedef typename function_space_t::elem_type_vector_t elem_type_vector_t;
	/** \brief integration result type */
	typedef Eigen::Matrix<typename kernel_t::scalar_t, Eigen::Dynamic, 1> result_vector_t;


protected:
	/**
	 * \brief evaluate integral on one element type and store the results into the result vector
	 * \details This is a helper functor called by tmp::call_each
	 * \tparam elem_t the element type over which integration is performed
	 */
	template <class elem_t>
	struct eval_on_elemtype
	{
		struct type
		{
			/** \brief quadrature order */
			static unsigned const num_quadrature_points = 5;
			/** \brief the field type */
			typedef field<elem_t, typename function_space_t::field_option> field_t;

			void operator() (weighted_integral_t &wi)
			{
				// initialise quadrature (this is done only once)
				gauss_quad<typename elem_t::domain_t, num_quadrature_points>::init();

				// integrate for each element of the same type
				std::for_each(
					wi.func_space.template begin<elem_t>(),
					wi.func_space.template end<elem_t>(),
					[&wi] (field_t const &f)
				{
					typedef typename weighted_field_integral<field_t, Kernel, num_quadrature_points>::result_t result_t;
					// autos are used to increase performance (reference is passed)
					result_t const &I = weighted_field_integral<field_t, Kernel, num_quadrature_points>::eval(f);
					typename field_t::dofs_t const &dofs = f.get_dofs();
					for (typename result_t::Index i = 0; i < field_t::num_dofs; ++i)
						wi.result_vector(dofs(i)) += I(i);
				}
				);
			}
		};
	};

public:
	/**
	 * \brief constructor initialises the function space reference member and allocates for the result
	 */
	weighted_integral(function_space_t const &func_space) : func_space(func_space)
	{
		result_vector.resize(func_space.get_num_dofs());
	}

	/**
	 * \brief evaluate integral and return reference to result vector
	 */
	result_vector_t const &eval(void)
	{
		result_vector = result_vector_t::Zero(result_vector.size());
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
};

#endif

