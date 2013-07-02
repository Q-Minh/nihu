/**
* \file function_space.hpp
* \ingroup funcspace
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class function_space
*/
#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "field.hpp"
#include "mesh.hpp"

/** \brief traits class of function spaces */
template <class Derived>
struct function_space_traits;


/** \brief CRTP base class of function spaces */
template <class Derived>
class function_space_base
{
private:
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

public:
	/** \brief the traits class */
	typedef function_space_traits<Derived> traits_t;

	/** \brief the field type vector */
	typedef typename traits_t::field_type_vector_t field_type_vector_t;

	/** \brief return begin iterator of a subvector of fields */
	template <class field_t>
	typename traits_t::template iterator<field_t>::type field_begin(void) const
	{
		return derived().template field_begin<field_t>();
	}

	/** \brief return end iterator of a subvector of fields */
	template <class field_t>
	typename traits_t::template iterator<field_t>::type field_end(void) const
	{
		return derived().template field_end<field_t>();
	}

	/** \brief return total number of degrees of freedoms */
	unsigned get_num_dofs(void) const
	{
		return derived().get_num_dofs();
	}
};


/** \brief forward declaration of function space view class */
template<class MeshT, class FieldOption>
class function_space_view;

/**
* \brief internal iterator class provides access to the mesh's elements as fields
* \tparam ElemType the element types that need to be accessed
*/
template <class FieldType>
class field_view_iterator_t : public mesh_elem_iterator_t<typename FieldType::elem_t>::type
{
	// CRTP check
	static_assert(std::is_base_of<field_base<FieldType>, FieldType>::value,
		"FieldType must be derived from field_base<FieldType>");
public:
	/** \brief the base iteartor type */
	typedef typename mesh_elem_iterator_t<typename FieldType::elem_t>::type base_it;
	typedef FieldType value_t;	/**< \brief the pointed data type */

	/**
	* \brief constructor from base iterator
	* \param it element iterator
	*/
	field_view_iterator_t(base_it const &it)
		: base_it(it)
	{
	}

	/**
	* \brief overloaded dereference operator simply converts dereferenced element into field
	* \return the referred field class
	*/
	value_t operator *(void) const
	{
		return value_t(base_it::operator*());
	}
};


/** \brief traits class of a function space view */
template <class MeshT, class FieldOption>
struct function_space_traits<function_space_view<MeshT, FieldOption> >
{
	/** \brief the field type vector */
	typedef typename tmp::transform<
		typename MeshT::elem_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1, tmp::_2> >,
		field_view<tmp::_1, FieldOption>
	>::type field_type_vector_t;

	/** \brief the iterator class traversing a field subvector */
	template <class field_t>
	struct iterator
	{
		typedef field_view_iterator_t<field_t> type;
	};
};


/**
* \brief helper class to return total number of degrees of freedom of a mesh
* \tparam mesh_t the mesh type
* \tparam field_option field generation option
*/
template <class mesh_t, class field_option>
struct get_num_dofs_impl;

/** \brief specialisation of get_num_dofs_impl for the constant case */
template <class mesh_t>
struct get_num_dofs_impl<mesh_t, constant_field>
{
	/**
	* \brief return number of degrees of freedom
	* \param mesh the mesh
	* \return total number of degrees of freedom
	*/
	static unsigned eval(mesh_t const &mesh)
	{
		return mesh.get_num_elements();
	}
};

/** \brief specialisation of get_num_dofs_impl for the isoparametric case */
template <class mesh_t>
struct get_num_dofs_impl<mesh_t, isoparametric_field>
{
	/**
	* \brief return number of degrees of freedom
	* \param mesh the mesh
	* \return total number of degrees of freedom
	*/
	static unsigned eval(mesh_t const &mesh)
	{
		return mesh.get_num_points();
	}
};



/**
* \brief FunctionSpace is a mesh extended with a Field generating option
* \tparam MeshT the underlying Mesh type
* \tparam FieldOption determines how the field is generated from the mesh
* \details The class is a proxy that stores a constant reference to the mesh.
* The class provides an iterator that can traverse the elements and derefers them as fields.
*/
template<class Mesh, class FieldOption>
class function_space_view : public function_space_base<function_space_view<Mesh, FieldOption> >
{
public:
	/** \brief the CRTP base class */
	typedef function_space_base<function_space_view<Mesh, FieldOption> > crtp_base;
	/** \brief the traits class */
	typedef typename crtp_base::traits_t traits_t;

	/** \brief the field type vector */
	typedef typename traits_t::field_type_vector_t field_type_vector_t;

	/** \brief template parameter as nested type */
	typedef Mesh mesh_t;
	/** \brief template parameter as nested type */
	typedef FieldOption field_option_t;

	/**
	* \brief constructor storing a reference from the mesh
	* \param mesh the mesh object to extend
	*/
	function_space_view(mesh_t const &mesh) : m_mesh(mesh)
	{
	}

	/**
	* \brief first field of given element type
	* \tparam ElemType the element type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type field_begin(void) const
	{
		// CRTP check
		static_assert(std::is_base_of<field_base<FieldType>, FieldType>::value,
			"FieldType must be derived from field_base<FieldType>");

		return m_mesh.template begin<typename FieldType::elem_t>(); // automatic conversion by iterator constructor
	}

	/**
	* \brief last field of given element type
	* \tparam ElemType the element type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type field_end(void) const
	{
		// CRTP check
		static_assert(std::is_base_of<field_base<FieldType>, FieldType>::value,
			"FieldType must be derived from field_base<FieldType>");

		return m_mesh.template end<typename FieldType::elem_t>(); // automatic conversion by iterator constructor
	}

	/**
	* \brief return number of degrees of freedom
	* \return number of degrees of freedom
	*/
	unsigned get_num_dofs(void) const
	{
		return get_num_dofs_impl<mesh_t, field_option_t>::eval(m_mesh);
	}

protected:
	/** \brief the stored mesh reference */
	mesh_t const &m_mesh;
};


/** \brief forward declaration of class function space */
template <class FieldTypeVector>
class function_space;

/** \brief Traits class of a function space */
template <class FieldTypeVector>
struct function_space_traits<function_space<FieldTypeVector> >
{
	/** \brief the field type vector */
	typedef FieldTypeVector field_type_vector_t;

	/** \brief the iterator type of a field type subvector */
	template <class field_t>
	struct iterator
	{
		typedef typename EigenStdVector<field_t>::type::const_iterator type;
	};
};


/** \brief metafunction to extract the element type of a field */
template <class field_t>
struct elemize
{
	typedef typename field_t::elem_t type;
};

/** \brief metafunction to return the element type vector of a field type vector */
template <class FieldTypeVector>
struct field_2_elem_type_vector
{
	typedef typename tmp::unique<typename tmp::transform<
		FieldTypeVector,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1, tmp::_2> >,
		elemize<tmp::_1>
	>::type>::type type;
};

/** \brief class describing a function space
* \tparam FieldTypeVector compile time vector of fields building the function space
*/
template <class FieldTypeVector>
class function_space :
	public function_space_base<function_space<FieldTypeVector> >,
	public mesh<typename field_2_elem_type_vector<FieldTypeVector>::type>
{
public:
	/** \brief the CRTP base class */
	typedef function_space_base<function_space<FieldTypeVector> > crtp_base;
	/** \brief the traits class */
	typedef typename crtp_base::traits_t traits_t;

	/** \brief the field type vector */
	typedef typename traits_t::field_type_vector_t field_type_vector_t;

	/** \brief the underlying mesh type */
	typedef mesh<typename field_2_elem_type_vector<FieldTypeVector>::type> mesh_t;
	
	/** \brief combine field_type_vector into a BIG heterogeneous std::vector container */
	typedef typename tmp::inherit<
		typename tmp::transform<
		field_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1,tmp::_2> >,
		EigenStdVector<tmp::_1>
		>::type
	>::type field_container_t;

protected:

	/** \brief subclass called by call_each to add a field to the function space */
	template <class field_t>
	struct field_adder { struct type {
		/** \brief add a field to the function space */
		bool operator() (unsigned const input[], function_space &fsp)
		{
			typedef typename field_t::elem_t elem_t;
			if (input[0] == field_t::id)
			{
				// construct element and push
				typename elem_t::nodes_t nodes;
				typename elem_t::coords_t coords;
				for (unsigned i = 0; i < elem_t::num_nodes; ++i)
				{
					nodes[i] = input[i+1];
					coords.col(i) = fsp.points[nodes[i]];
				}
				elem_t const &elemref = fsp.push_element(elem_t(coords, fsp.m_num_elements++, nodes));

				// construct field and push
				typename field_t::dofs_t dofs;
				for (unsigned i = 0; i < field_t::num_dofs; ++i)
					dofs[i] = input[i+elem_t::num_nodes+1];
				fsp.push_field(field_t(elemref, dofs));

				fsp.m_num_dofs = std::max(fsp.m_num_dofs, dofs.maxCoeff()+1);
                return true;
			}
			return false;
		}
	};};

private:
	/** \brief disable copying of function spaces as the fields contain references to the elements */
	function_space(function_space const &other);

	/** \brief disable assignment of function spaces as the fields contain references to the elements */
	function_space const &operator=(function_space const &other);


public:
	/** \brief constructor */
	function_space() : m_num_dofs(0)
	{
	}

	/** \brief constructor from node and field definition matrices */
	template <class node_matrix_t, class field_matrix_t>
	function_space(node_matrix_t const &nodes, field_matrix_t const &fields)
		: m_num_dofs(0)
	{
		unsigned const N_MAX_FIELD = 1000;
		double c[mesh_t::nDim];

		for (int i = 0; i < nodes.rows(); ++i)
		{
			for (unsigned j = 0; j < mesh_t::nDim; ++j)
				c[j] = nodes(i,j);
			this->add_node(c);
		}
		unsigned e[N_MAX_FIELD];
		for (int i = 0; i < fields.rows(); ++i)
		{
			for (int j = 0; j < fields.cols(); ++j)
				e[j] = (unsigned)fields(i,j);
			add_field(e);
		}
	}

	/** \brief return begin iterator of a subspace */
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_begin(void) const
	{
		return m_fields.EigenStdVector<FieldType>::type::begin();
	}

	/** \brief return end iterator of a subspace */
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_end(void) const
	{
		return m_fields.EigenStdVector<FieldType>::type::end();
	}

	/** \brief add a field to the function space */
	bool add_field(unsigned const input[])
	{
		return tmp::call_until<
			field_type_vector_t,
			field_adder<tmp::_1>,
			unsigned const*,
			function_space &
		>(input, *this);
	}

	/** \brief push a field to the vector of fields */
	template <class field_t>
	field_t const &push_field(field_t const &f)
	{
		m_fields.EigenStdVector<field_t>::type::push_back(f);
		return *(m_fields.EigenStdVector<field_t>::type::rbegin());
	}

	/**
	 * \brief return number of fields
	 * \return number of fields in the function space
	 */
	unsigned get_num_fields(void) const
	{
		return this->get_num_elements();
	}

	/** \brief return number of dofs */
	unsigned get_num_dofs(void) const
	{
		return m_num_dofs;
	}


protected:
	/** \brief fields (BIG heterogeneous container) */
	field_container_t m_fields;	
	/** \brief number of degrees of freedoms */
	unsigned m_num_dofs;
};


#endif

