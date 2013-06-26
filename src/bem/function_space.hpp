/**
* \file function_space.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class function_space
*/
#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "field.hpp"
#include "mesh.hpp"

template <class Derived>
struct function_space_traits;


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
	typedef function_space_traits<Derived> traits_t;

	template <class field_t>
	typename traits_t::template iterator<field_t>::type field_begin(void) const
	{
		return derived().template field_begin<field_t>();
	}

	template <class field_t>
	typename traits_t::template iterator<field_t>::type field_end(void) const
	{
		return derived().template field_end<field_t>();
	}
};


template<class MeshT, class FieldOption>
class function_space_view;

/**
* \brief internal iterator class provides access to the mesh's elements as fields
* \tparam ElemType the element types that need to be accessed
* \todo this iterator should only be used in field view context
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


template <class MeshT, class FieldOption>
struct function_space_traits<function_space_view<MeshT, FieldOption> >
{
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
template<class MeshT, class FieldOption>
class function_space_view : public function_space_base<function_space_view<MeshT, FieldOption> >
{
public:
	/** \brief template parameter as nested type */
	typedef MeshT mesh_t;
	/** \brief template parameter as nested type */
	typedef FieldOption field_option_t;

	/** \brief elem_type_vector inherited from mesh */
	typedef typename mesh_t::elem_type_vector_t elem_type_vector_t;

	/** \brief metafunction to convert an element type into a field type */
	template <class elem_t>
	struct fieldize
	{
		typedef field_view<elem_t, field_option_t> type; /**< \brief metafunction return type */
	};

	/** \brief a vector of field types computed from the element type vector */
	typedef typename tmp::transform<
		elem_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1, tmp::_2> >,
		fieldize<tmp::_1>
	>::type field_type_vector_t;

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
	field_view_iterator_t<FieldType> field_begin(void) const
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
	field_view_iterator_t<FieldType> field_end(void) const
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



















template <class FieldType>
struct function_space_field_iterator_t
{
	typedef typename EIGENSTDVECTOR(FieldType)::const_iterator type;
};


template <class FieldTypeVector>
class function_space;



template <class FieldTypeVector>
struct function_space_traits<function_space<FieldTypeVector> >
{
	template <class field_t>
	struct iterator
	{
		typedef function_space_field_iterator_t<field_t> type;
	};
};


/** \brief metafunction computing the first field's x_t in a vector of fields */
template <class FieldTypeVector>
struct first_fields_x_type
{
	typedef typename tmp::deref<
		typename tmp::begin<FieldTypeVector>::type
	>::type::elem_t::x_t type;
};

/**
 * \brief container class for a function_space
 * \tparam FieldTypeVector compile time vector of the contained field types
 */
template <class FieldTypeVector>
class function_space :
	public function_space_base<function_space<FieldTypeVector> >,
	public field_points<typename first_elements_x_type<FieldTypeVector>::type>
{
public:
	/** \brief define template parameter as nested type */
	typedef FieldTypeVector field_type_vector_t;

	/** \brief type of base class */
	typedef field_points<typename first_fields_x_type<FieldTypeVector>::type> base_t;

	static unsigned const nDim = base_t::nDim;	/**< \brief number of dimensions of the mesh */

	/** \brief metafunction to convert T into std::vector<T> */
	template <class T>
	struct vectorize { typedef EIGENSTDVECTOR(T) type; };

	/** \brief combine field_vector into a BIG heterogeneous std::vector container */
	typedef typename tmp::inherit<
		typename tmp::transform<
		field_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1,tmp::_2> >,
		vectorize<tmp::_1>
		>::type
	>::type field_container_t;

	/** \brief type of a nodal vector */
	typedef typename base_t::x_t x_t;

	/** \brief type of a nodal coordinate */
	typedef typename x_t::Scalar scalar_t;

	template <class FieldType>
	struct field_iterator_t
	{
		typedef typename function_space_field_iterator_t<FieldType>::type type;
	};


protected:
	template <class field_t>
	struct field_adder { struct type	{
		/**
		 * \brief add a field of given type to the function space
		 * \param m the function space to extend
		 * \param input array containing element node indices
		 */
		bool operator() (unsigned const input[], function_space &fs)
		{
			if (input[0] == field_t::field_id)
			{
				// construct element
				typename field_t::elem_t::nodes_t nodes;
				typename field_t::elem_t::dofs_t dofs;
				typename field_t::elem_t::coords_t coords;
				for (unsigned i = 0; i < field_t::elem_t::num_nodes; ++i)
				{
					nodes[i] = input[i+1];
					coords.row(i) = fs.points[nodes[i]];
				}
				for (unsigned i = 0; i < field_t::nset_t::num_nodes; ++i)
					dofs[i] = input[i+field_t::elem_t::num_nodes+1];
				fs.push_field(field_t(elem_t(coords, fs.m_num_fields++, nodes), dofs));
                return true;
			}
			return false;
		}
	};};

public:
	function_space() : m_num_fields(0)
	{
	}

	/**
	 * \brief build the mesh from MATLAB matrices
	 * \tparam type of nodes matrix
	 * \tparam type of elements matrix
	 * \param [in] nodes matrix of nodal coordinates
	 * \param [in] elements matrix of element node indices
	 */
	template <class node_t, class field_t>
	function_space(node_t const &nodes, field_t const &elements) : m_num_elements(0)
	{
		unsigned const N_MAX_ELEM = 100;
		double c[nDim];

		for (int i = 0; i < nodes.rows(); ++i)
		{
			for (unsigned j = 0; j < nDim; ++j)
				c[j] = nodes(i,j);
			add_node(c);
		}
		unsigned e[N_MAX_ELEM];
		for (int i = 0; i < elements.rows(); ++i)
		{
			for (int j = 0; j < elements.cols(); ++j)
				e[j] = (unsigned)elements(i,j);
			add_field(e);
		}
	}

	/**
	 * \brief return begin iterator of the elements
	 */
	template <class FieldType>
	typename field_iterator_t<FieldType>::type begin(void) const
	{
		return m_fields.EIGENSTDVECTOR(FieldType)::begin();
	}

	/**
	 * \brief return end iterator of the elements
	 */
	template <class FieldType>
	typename field_iterator_t<FieldType>::type end(void) const
	{
		return m_fields.EIGENSTDVECTOR(ElemType)::end();
	}

	/**
	 * \brief add a new element to the mesh
	 * \param input array of unsigned values.
	 * input[0] is the elem ID, the subsequent elements are the nodal indices in the mesh
	 * \return true if the element is inserted into the mesh
	 */
	bool add_field(unsigned const input[])
	{
		return tmp::call_until<
			field_type_vector_t,
			field_adder<tmp::_1>,
			unsigned const*,
			function_space &
		>(input, *this);
	}

	/**
	 * \brief add a new node to the mesh
	 * \param input array of scalars containing the coordinates
	 */
	void add_node(scalar_t input[])
	{
		x_t c;
		for (unsigned i = 0; i < nDim; ++i)
			c[i] = input[i];
		this->add_point(c);
	}

	/**
	 * \brief add a new element to the mesh
	 * \tparam elem_t the element type
	 * \param e the element to be added
	 */
	template <class field_t>
	void push_field(field_t const &e)
	{
		m_fields.EIGENSTDVECTOR(field_t)::push_back(e);
	}

	/**
	 * \brief return number of elements
	 * \return number of elements in the mesh
	 */
	unsigned get_num_fields(void) const
	{
		return m_num_fields;
	}

protected:
	field_container_t m_fields;	/**< \brief element geometries (BIG heterogeneous container) */
	unsigned m_num_fields;	/**< \brief total number of elements in the mesh */
};



#endif

