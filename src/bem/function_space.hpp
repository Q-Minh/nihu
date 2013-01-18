/**
 * \file function_space.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class FunctionSpace
 */
#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "field.hpp"
#include "mesh.hpp"

template<class ElemVector, class FieldOption>
class FunctionSpace
{
public:
	typedef ElemVector elem_type_vector_t;
	typedef FieldOption field_option;
	
	typedef Mesh<elem_type_vector_t> mesh_t;

	FunctionSpace(mesh_t const &mesh) : mesh(mesh)
	{
	}
	
	void go(void)
	{
		std::for_each(
			mesh.template elembegin<tria_1_elem>(),
			mesh.template elemend<tria_1_elem>(),
			[] (tria_1_elem const &e)
			{
				tria_1_elem::xi_t xi = tria_1_elem::xi_t::Zero();
				double f = e.get_x(xi).norm();
				typename Field<tria_1_elem, field_option>::nset_t::L_t N = Field<tria_1_elem, field_option>::nset_t::eval_L(xi);
				std::cout << f * N << std::endl;
			}
		);
	}

protected:
	mesh_t const &mesh;
};

#endif

