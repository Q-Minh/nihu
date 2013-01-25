/**
 * \file rayleigh.hpp
 * \brief implementation of the Rayleigh integeral in 3D
 */
 #ifndef BEM_HPP_INCLUDED
 #define BEM_HPP_INCLUDED

 #include "weighted_integral.hpp"

/**
 * \brief evaluate the BEM radiation integral on a function space
 * \details The integral is evaluated over a function space
 * consisting of elements listed in the template parameter ElemVector and
 * extended into a function space defined by FieldOption.
 * \tparam ElemVector vector of element types that can be contained by the mesh
 * \tparam FieldOption function space generation option (constant_field or isoparametric_field)
 */
template <class ElemVector, class FieldOption>
class bem
{
public:
	/** \brief template parameter as nested type */
	typedef ElemVector elem_type_vector_t;
	/** \brief template parameter as nested type */
	typedef FieldOption field_option_t;

	typedef Mesh<elem_type_vector_t> mesh_t;
	typedef typename mesh_t::x_t x_t;
	typedef function_space<mesh_t, field_option_t> function_space_t;
	typedef green_HG_kernel kernel_t;
	typedef weighted_integral<function_space_t, kernel_t> weighted_integral_t;
	typedef typename weighted_integral_t::result_vector_t result_vector_t;

	/**
	 * \brief the constructor initialises the stored function space and weighted integral
	 * \param mesh the stored mesh
	 */
	bem(mesh_t const &mesh) : f_space(mesh), wi(f_space)
	{
	}

	/**
	 * \brief evaluate the bem integral for a source point and a wave number
	 * \param x0 the source point
	 * \param the k wave number
	 */
	result_vector_t const &eval(x_t const &x0, dcomplex const &k)
	{
		kernel_t::set_x0(x0);
		kernel_t::set_wave_number(k);
		return wi.eval();
	}

protected:
	/** \brief the function space object over which integration is performed */
	function_space_t f_space;
	/** \brief the weighted integral object that performs integration */
	weighted_integral_t wi;
};

#endif

