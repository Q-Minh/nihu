/**
 * \file rayleigh.hpp
 * \brief implementation of the Rayleigh integeral in 3D
 */
 #ifndef RAYLEIGH_HPP_INCLUDED
 #define RAYLEIGH_HPP_INCLUDED

 #include "weighted_integral.hpp"

template <class ElemVector, class FieldOption>
class rayleigh
{
public:
	/** \brief template parameter as nested type */
	typedef ElemVector elem_type_vector_t;
	/** \brief template parameter as nested type */
	typedef FieldOption field_option_t;

	typedef Mesh<elem_type_vector_t> mesh_t;
	typedef typename mesh_t::x_t x_t;
	typedef function_space<mesh_t, field_option_t> function_space_t;
	typedef green_kernel kernel_t;
	typedef weighted_integral<function_space_t, kernel_t> weighted_integral_t;
	typedef typename weighted_integral_t::result_vector_t result_vector_t;

	rayleigh(mesh_t const &mesh) : f_space(mesh), wi(f_space)
	{
	}

	result_vector_t const &eval(x_t const &x0, dcomplex const &k)
	{
		kernel_t::set_x0(x0);
		kernel_t::set_wave_number(k);
		return wi.eval();
	}

protected:
	function_space_t f_space;
	weighted_integral_t wi;
};

#endif

