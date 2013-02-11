/**
 * \file bem.hpp
 * \brief implementation of the Helmholtz BEM radiation problem in 3D
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
template <class ElemVector, class FieldOption, class Kernel>
class bem
{
public:
	/** \brief template parameter as nested type */
	typedef ElemVector elem_type_vector_t;
	/** \brief template parameter as nested type */
	typedef FieldOption field_option_t;
	/** \brief specify the actual kernel to integrate */
	typedef Kernel kernel_t;

	/** \brief the stored mesh type */
	typedef Mesh<elem_type_vector_t> mesh_t;
	/** \brief the location vector type */
	typedef typename mesh_t::x_t x_t;
	/** \brief the function space type */
	typedef function_space<mesh_t, field_option_t> function_space_t;
	/** \brief the weighted integral type */
	typedef weighted_integral<function_space_t, kernel_t> weighted_integral_t;
	/** \brief the result type of integration */
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
	 * \param k the wave number
	 * \return reference to the internally stored result vector
	 */
	result_vector_t const &eval(x_t const &x0, dcomplex const &k)
	{
		// initialise the kernel
		kernel_t::set_x0(x0);
		kernel_t::set_wave_number(k);
		// evaluate the weighted integral
		return wi.eval();
	}

protected:
	/** \brief the function space object over which integration is performed */
	function_space_t f_space;
	/** \brief the weighted integral object that integrates and stores the result */
	weighted_integral_t wi;
};

#endif

