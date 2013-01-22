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

	typedef typename first_elements_x_type<elem_type_vector_t>::type x_t;
	typedef field_points<x_t> field_points_t;
	typedef Mesh<elem_type_vector_t> mesh_t;
	typedef function_space<mesh_t, field_option_t> function_space_t;
	typedef green_kernel kernel_t;

	rayleigh(mesh_t const &mesh, field_points_t const &field_p, dcomplex wave_number)
	: mesh(mesh), field_p(field_p), wave_number(wave_number)
	{
	}

	void eval(void)
	{
		kernel_t::set_wave_number(wave_number);

		std::for_each(
			field_p.begin(),
			field_p.end(),
			[this] (x_t const &x) {
				kernel_t::set_x0(x);
				weighted_integral::eval<function_space_t, kernel_t>(function_space_t(mesh));
			}
		);
	}

protected:
	mesh_t mesh;
	field_points_t field_p;
	dcomplex wave_number;
};

#endif
