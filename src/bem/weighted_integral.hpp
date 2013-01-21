/**
* \file weighted_integral.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class weighted_integral
*/
#ifndef WEIGHTED_INTEGRAL_SPACE_HPP_INCLUDED
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
			xi_t xi = xivec.row(i);
			kernel_input_t input(field.get_elem(), xi);
			kernel_result_t kernel_res = kernel_t::eval(input);
			I += nset_t::eval_L(xi) * (kernel_res * input.get_jacobian() * weightvec(i));
		}

		return I;
	}
};


template<class FunctionSpace, class Kernel>
class weighted_integral
{
public:
	typedef FunctionSpace function_space_t;
	typedef Kernel kernel_t;

	typedef typename function_space_t::mesh_t mesh_t;
	typedef typename mesh_t::elem_type_vector_t elem_type_vector_t;

	weighted_integral(function_space_t const &func_space) : func_space(func_space)
	{
	}

	void do_it(void)
	{
		typedef tria_1_elem elem_t;
		typedef field<elem_t, typename function_space_t::field_option> field_t;
		std::for_each(
			func_space.get_mesh().template elembegin<elem_t>(),
			func_space.get_mesh().template elemend<elem_t>(),
			[] (elem_t const &e)
		{
			std::cout << weighted_field_integral<field_t, kernel_t, 5>::eval(field_t(e)) << std::endl;
		}
		);
	}

protected:
	function_space_t const &func_space;
};

#endif

