/**
 * \file function_space.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class function_space
 */
#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "field.hpp"
#include "mesh.hpp"

/**
 * \brief FunctionSpace is a mesh extended with a Field generating option
 * \tparam MeshT the underlying Mesh type
 * \tparam FieldOption determines how the field is generated from the mesh
 * \details The class is a proxy that stores a constant reference to the mesh.
 * The class provides an iterator that can traverse the elements and derefers them as fields.
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
	 * \brief internal iterator class provides access to the mesh's elements as fields
	 * \tparam ElemType the element types that need to be accessed
	 */
	template <class ElemType>
	class field_iterator_t : public mesh_t::template elem_iterator_t<ElemType>::type
	{
	public:
		typedef typename mesh_t::template elem_iterator_t<ElemType>::type base_it;
		typedef field<ElemType, field_option> value_t;

		/**
		 * \brief constructor from base iterator
		 * \param it element iterator
		 */
		field_iterator_t(base_it const &it) : base_it(it) { }

		/**
		 * \brief overloaded dereference operator simply converts dereferenced element into field
		 * \return the referred field class
		 */
		value_t operator *(void) const
		{
			return value_t(base_it::operator*());
		}
	};

	/**
	 * \brief constructor storing a reference from the mesh
	 * \param mesh the mesh object to extend
	 */
	function_space(mesh_t const &mesh) : mesh(mesh)
	{
	}

	/**
	 * \brief first field of given element type
	 * \tparam ElemType the element type to access
	 */
	template <class ElemType>
	field_iterator_t<ElemType> begin(void) const
	{
		return mesh.begin<ElemType>(); // automatic conversion by iterator constructor
	}

	/**
	 * \brief last field of given element type
	 * \tparam ElemType the element type to access
	 */
	template <class ElemType>
	field_iterator_t<ElemType> end(void) const
	{
		return mesh.end<ElemType>(); // automatic conversion by iterator constructor
	}

protected:
	/** \brief the stored mesh reference */
	mesh_t const &mesh;
};

#endif

