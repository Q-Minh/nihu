/**
* \file function_space.hpp
* \ingroup funcspace
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class function_space
*/
#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "../util/casted_iterator.hpp"
#include "field.hpp"
#include "mesh.hpp"

/** \brief traits class of function spaces */
template <class Derived>
struct function_space_traits;


/** \brief CRTP base class of function spaces */
template <class Derived>
class function_space_base
{
public:
	NIHU_CRTP_HELPERS

	/** \brief the traits class */
	typedef function_space_traits<Derived> traits_t;

	/** \brief the field type vector */
	typedef typename traits_t::field_type_vector_t field_type_vector_t;
};


template <class Derived>
class function_space_impl;


// forward declaration
template<class Mesh, class FieldOption>
class function_space_view;


/** \brief traits class of a function space view */
template <class Mesh, class FieldOption>
struct function_space_traits<function_space_view<Mesh, FieldOption> >
{
	/** \brief the field type vector */
	typedef typename tmp::transform<
		typename Mesh::elem_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1, tmp::_2> >,
		field_view<tmp::_1, FieldOption>
	>::type field_type_vector_t;

	/** \brief the iterator class traversing a field subvector */
	template <class field_t>
	struct iterator
	{
		typedef casted_iterator<
			typename mesh_elem_iterator_t<typename field_t::elem_t>::type,
			field_t
		> type;
	};
};


/**
* \brief function_space_view is a mesh extended with a Field generating option
* \tparam Mesh the underlying Mesh type
* \tparam FieldOption determines how the field is generated from the mesh
* The class provides an iterator that can traverse the elements and derefers them as fields.
*/
template<class Mesh, class FieldOption>
class function_space_impl<function_space_view<Mesh, FieldOption> > :
	public Mesh
{
public:
	/** \brief the traits class */
	typedef function_space_traits< function_space_view<Mesh, FieldOption> > traits_t;

	/** \brief template parameter as nested type */
	typedef Mesh mesh_t;
	/** \brief template parameter as nested type */
	typedef FieldOption field_option_t;

	/**
	* \brief first field of given element type
	* \tparam ElemType the element type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_begin(void) const
	{
		return mesh_t::template begin<typename FieldType::elem_t>();
	}

	/**
	* \brief last field of given element type
	* \tparam ElemType the element type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_end(void) const
	{
		return mesh_t::template end<typename FieldType::elem_t>();
	}

private:
	unsigned get_num_dofs_impl(field_option::constant) const
	{
		return mesh_t::get_num_elements();
	}

	unsigned get_num_dofs_impl(field_option::isoparametric) const
	{
		return mesh_t::get_num_points();
	}


public:
	/**
	* \brief return number of degrees of freedom
	* \return number of degrees of freedom
	*/
	unsigned get_num_dofs(void) const
	{
		return get_num_dofs_impl(field_option_t());
	}
};



/**
* \brief function_space_view is a mesh extended with a Field generating option
* \tparam Mesh the underlying Mesh type
* \tparam FieldOption determines how the field is generated from the mesh
* The class provides an iterator that can traverse the elements and derefers them as fields.
*/
template<class Mesh, class FieldOption>
class function_space_view :
	public function_space_base<function_space_view<Mesh, FieldOption> >,
	public function_space_impl<function_space_view<Mesh, FieldOption> >
{
};


template <class Mesh, class Option>
function_space_view<Mesh, Option> const &
	create_function_space_view(Mesh const &msh, Option)
{
	return static_cast<function_space_view<Mesh, Option> const &>(msh);
}

template <class Mesh>
function_space_view<Mesh, field_option::isoparametric> const &
	isoparametric_view(Mesh const &msh)
{
	return create_function_space_view(msh, field_option::isoparametric());
}

template <class Mesh>
function_space_view<Mesh, field_option::constant> const &
	constant_view(Mesh const &msh)
{
	return create_function_space_view(msh, field_option::constant());
}




template <class FuncSpace>
class dirac_space;

template <class FuncSpace>
struct function_space_traits<dirac_space<FuncSpace> >
{
	/** \brief the field type vector */
	typedef typename tmp::transform<
		typename function_space_traits<FuncSpace>::field_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1, tmp::_2> >,
		dirac_field<tmp::_1>
	>::type field_type_vector_t;

	/** \brief the iterator class traversing a field subvector */
	template <class dirac_field_t>
	struct iterator
	{
		typedef casted_iterator<
			typename function_space_traits<FuncSpace>::template iterator<
				typename dirac_field_t::field_t
			>::type,
			dirac_field_t,
			field_impl<typename dirac_field_t::field_t>
		> type;
	};
};


template <class FuncSpace>
class dirac_space :
	public function_space_base<dirac_space<FuncSpace> >,
	public function_space_impl<FuncSpace>
{
public:
	/** \brief the traits class */
	typedef function_space_traits<dirac_space<FuncSpace> > traits_t;
	typedef function_space_impl<FuncSpace> impl_t;

	/** \brief the field type vector */
	typedef typename traits_t::field_type_vector_t field_type_vector_t;

	/** \brief return begin iterator of a subvector of fields */
	template <class dirac_field_t>
	typename traits_t::template iterator<dirac_field_t>::type
		field_begin(void) const
	{
		return impl_t::template field_begin<typename dirac_field_t::field_t>();
	}

	/** \brief return end iterator of a subvector of fields */
	template <class dirac_field_t>
	typename traits_t::template iterator<dirac_field_t>::type
		field_end(void) const
	{
		return impl_t::template field_end<typename dirac_field_t::field_t>();
	}
};


/** \brief factory function to convert a function space into a dirac space */
template <class FuncSpace>
dirac_space<FuncSpace> const &
	dirac(function_space_base<FuncSpace> const &space)
{
	return static_cast<dirac_space<FuncSpace> const &>(
		static_cast<function_space_impl<FuncSpace> const &>(space.derived()));
}



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
class function_space_impl<function_space<FieldTypeVector> > :
	public mesh<typename field_2_elem_type_vector<FieldTypeVector>::type>
{
public:
	/** \brief the traits class */
	typedef typename function_space_traits<function_space<FieldTypeVector> >::traits_t traits_t;

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
		bool operator() (unsigned const input[], function_space_impl &fsp)
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
	function_space_impl(function_space_impl const &other);

	/** \brief disable assignment of function spaces as the fields contain references to the elements */
	function_space_impl const &operator=(function_space_impl const &other);


public:
	/** \brief constructor */
	function_space_impl() : m_num_dofs(0)
	{
	}

	/** \brief constructor from node and field definition matrices */
	template <class node_matrix_t, class field_matrix_t>
	function_space_impl(node_matrix_t const &nodes, field_matrix_t const &fields)
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
			function_space_impl &
		>(input, *this);
	}

	/** \brief push a field to the vector of fields */
	template <class field_t>
	field_t const &push_field(field_base<field_t> const &f)
	{
		m_fields.EigenStdVector<field_t>::type::push_back(f.derived());
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








/** \brief class describing a function space
* \tparam FieldTypeVector compile time vector of fields building the function space
*/
template <class FieldTypeVector>
class function_space :
	public function_space_base<function_space<FieldTypeVector> >,
	public function_space_impl<function_space<FieldTypeVector> >
{
public:
	typedef function_space_impl<function_space<FieldTypeVector> > impl_t;

	/** \brief constructor */
	function_space() : impl_t()
	{
	}

	/** \brief constructor from node and field definition matrices */
	template <class node_matrix_t, class field_matrix_t>
	function_space(node_matrix_t const &nodes, field_matrix_t const &fields)
		: impl_t(nodes, fields)
	{
	}
};


#endif

