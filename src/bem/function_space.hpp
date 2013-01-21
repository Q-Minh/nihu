/**
 * \file function_space.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class FunctionSpace
 */
#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "field.hpp"
#include "mesh.hpp"

/**
 * \brief FunctionSpace is a mesh extended with a Field generating option
 * \tparam MeshT the underlying Mesh type
 * \tparam FieldOption determines how the field is generated from the mesh
 */
template<class MeshT, class FieldOption>
class function_space
{
public:
	/** \brief template parameter as nested type */
	typedef MeshT mesh_t;
	/** \brief template parameter as nested type */
	typedef FieldOption field_option;
	
	/** \brief elem_type_vector inherited from mesh */
	typedef typename mesh_t::elem_type_vector_t elem_type_vector_t;

	/**
	 * \brief constructor storing a reference from the mesh
	 * \param mesh the mesh object to extend
	 */
	function_space(mesh_t const &mesh) : mesh(mesh)
	{
	}

	mesh_t const &get_mesh(void) const
	{
		return mesh;
	}
	
	void go(void) const
	{
		std::for_each(
			mesh.template elembegin<tria_1_elem>(),
			mesh.template elemend<tria_1_elem>(),
			[] (tria_1_elem const &e)
			{
				tria_1_elem::xi_t xi = tria_1_elem::xi_t::Zero();
				double f = e.get_x(xi).norm();
				typedef typename field<tria_1_elem, field_option>::nset_t nset_t;
				typedef typename nset_t::L_t L_t;
				L_t N = nset_t::eval_L(xi);
				std::cout << f * N << std::endl;
			}
		);
	}

protected:
	/** \brief the stored mesh reference */
	mesh_t const &mesh;
};

#endif

