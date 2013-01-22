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

template <class Field, class Kernel, unsigned N>
class weighted_field_integral
{
public:
	typedef Field field_t;
	typedef Kernel kernel_t;

	typedef typename kernel_t::input_t kernel_input_t;
	typedef typename kernel_t::result_t kernel_result_t;

	typedef typename field_t::elem_t elem_t;
	typedef typename elem_t::xi_t xi_t;
	typedef typename elem_t::domain_t domain_t;
	typedef gauss_quad<domain_t, N> quadrature_t;
	typedef typename quadrature_t::xivec_t xivec_t;
	typedef typename quadrature_t::weightvec_t weightvec_t;
	typedef typename field_t::nset_t nset_t;
	typedef typename Eigen::Matrix<dcomplex, field_t::nset_t::num_nodes, 1> result_t;

	static result_t eval(field_t const &field)
	{
		xivec_t xivec = quadrature_t::get_xi();
		weightvec_t weightvec = quadrature_t::get_weight();

		result_t I = result_t();
		for (typename xivec_t::Index i = 0; i < quadrature_t::size; ++i)
		{
			auto xi = xivec.row(i);
			kernel_input_t input(field.get_elem(), xi);
			kernel_result_t kernel_res = kernel_t::eval(input);
			I += nset_t::eval_L(xi) * (kernel_res * input.get_jacobian() * weightvec(i));
		}

		return I;
	}
};


/**
 * \brief 
 */
template <class FunctionSpace, class Kernel>
class weighted_integral
{
	typedef typename FunctionSpace::elem_type_vector_t elem_type_vector_t;
		
	template <class elem_t>
	struct eval_on_elemtype
	{
		struct type
		{
			void operator() (FunctionSpace const &func_space)
			{
				gauss_quad<typename elem_t::domain_t, 5>::init();
				typedef field<elem_t, typename FunctionSpace::field_option> field_t;
				
				std::for_each(
					func_space.template begin<elem_t>(),
					func_space.template end<elem_t>(),
					[] (field_t const &f)
					{
						std::cout << weighted_field_integral<field_t, Kernel, 5>::eval(f) << std::endl << std::endl;
					}
				);
			}
		};
	};

public:
	static void eval(FunctionSpace const &func_space)
	{
		tmp::call_each<
			elem_type_vector_t,
			eval_on_elemtype<tmp::_1>,
			FunctionSpace const&
		>(func_space);
	}
};

#endif

