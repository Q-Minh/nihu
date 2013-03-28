/**
 * \file bem.hpp
 * \brief implementation of the radiation problem in 3D
 */
#ifndef BEM_HPP_INCLUDED
#define BEM_HPP_INCLUDED

#include "weighted_integral.hpp"
#include "weighted_surface_integral.hpp"

/**
 * \brief evaluate the BEM radiation integral on a function space
 * \details The integral is evaluated over a function space
 * consisting of elements listed in the template parameter ElemVector and
 * extended into a function space defined by FieldOption.
 * \tparam ElemVector vector of element types that can be contained by the mesh
 * \tparam FieldOption function space generation option (constant_field or isoparametric_field)
 */
template <class Kernel, class ElemVector, class FieldOption>
class bem_radiation
{
public:
	typedef Kernel kernel_t; ///< template argument as nested type
	typedef ElemVector elem_type_vector_t; 	///< \brief template parameter as nested type
	typedef FieldOption field_option_t; 	///< \brief template parameter as nested type


	typedef Mesh<elem_type_vector_t> mesh_t;	///< \brief the stored mesh type
	typedef typename mesh_t::x_t x_t;			///< \brief the location vector type
	typedef function_space<mesh_t, field_option_t> function_space_t;	///< \brief the function space type
	typedef weighted_integral<kernel_t, function_space_t> weighted_integral_t;	///< \brief the weighted integral type
	typedef typename weighted_integral_t::result_vector_t result_vector_t;	///< \brief the result type of integration

	/**
	 * \brief the constructor initialises the stored function space and weighted integral
	 * \param mesh the stored mesh
	 */
	bem_radiation(mesh_t const &mesh) : f_space(mesh), wi(f_space) {}

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
		return wi.integrate();
	}

protected:
	function_space_t f_space;	/// < \brief the function space over which integration is performed
	weighted_integral_t wi;		///< \brief object that integrates and stores the result
};


/**
 * \brief evaluate the BEM radiation integral on a function space
 * \details The integral is evaluated over a function space
 * consisting of elements listed in the template parameter ElemVector and
 * extended into a function space defined by FieldOption.
 * \tparam ElemVector vector of element types that can be contained by the mesh
 * \tparam FieldOption function space generation option (constant_field or isoparametric_field)
 */
template <class Kernel, class ElemVector, class FieldOption1, class FieldOption2>
class bem_surface_system
{
public:
	typedef Kernel kernel_t;	/** \brief template parameter as nested type */
	typedef ElemVector elem_type_vector_t;	/** \brief template parameter as nested type */
	typedef FieldOption1 field_option1_t;	/** \brief template parameter as nested type */
	typedef FieldOption2 field_option2_t;	/** \brief template parameter as nested type */

	typedef Mesh<elem_type_vector_t> mesh_t;
	typedef function_space<mesh_t, field_option1_t> function_space1_t;
	typedef function_space<mesh_t, field_option2_t> function_space2_t;
	typedef weighted_surface_integral<kernel_t, function_space1_t, function_space2_t> weighted_surface_integral_t;

	/**
	 * \brief the constructor initialises the stored function space and weighted integral
	 * \param mesh the stored mesh
	 */
	bem_surface_system(mesh_t const &mesh)
		: m_f_space1(mesh), m_f_space2(mesh), m_wsi(m_f_space1, m_f_space2)
	{
	}

	/**
	 * \brief evaluate the bem integral for a source point and a wave number
	 * \param x0 the source point
	 * \param k the wave number
	 * \return reference to the internally stored result vector
	 */
	template <class result_matrix_t>
	void eval(dcomplex const &k, result_matrix_t &result_matrix)
	{
		// initialise the kernel
		kernel_t::set_wave_number(k);
		// evaluate the weighted integral
		m_wsi.integrate(result_matrix);
	}

protected:
	function_space1_t m_f_space1;	///< \brief function space over which integration is performed
	function_space2_t m_f_space2;	///< \brief function space over which integration is performed

	weighted_surface_integral_t m_wsi;	///< \brief the object that integrates
};

#endif
