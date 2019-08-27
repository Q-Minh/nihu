// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
* \file function_space.hpp
* \ingroup funcspace
* \brief declaration of class function_space
*/
#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "../util/casted_iterator.hpp"
#include "field.hpp"
#include "mesh.hpp"

namespace NiHu
{

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

	/** \brief the underlying mesh type */
	typedef typename traits_t::mesh_t mesh_t;

	/** \brief the field type vector */
	typedef typename traits_t::field_type_vector_t field_type_vector_t;

	/** \brief return reference to the underlying mesh
	* \return reference to the mesh part */
	mesh_t const &get_mesh() const
	{
		return derived().get_mesh();
	}

	/**
	* \brief begin iterator of given field type
	* \tparam FieldType the field type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_begin() const
	{
		return derived().template field_begin<FieldType>();
	}

	/**
	* \brief end iterator of given field type
	* \tparam FieldType the field type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_end() const
	{
		return derived().template field_end<FieldType>();
	}

	/**
	* \brief return number of degrees of freedom
	* \return number of degrees of freedom
	*/
	size_t get_num_dofs() const
	{
		return derived().get_num_dofs();
	}
};


/** \brief implementation class of function spaces */
template <class Derived>
class function_space_impl;


// forward declaration
template<class Mesh, class FieldOption, class Dimension>
class function_space_view;


/** \brief traits class of a function space view */
template <class Mesh, class FieldOption, class Dimension>
struct function_space_traits<function_space_view<Mesh, FieldOption, Dimension> >
{
	/** \brief the underlying mesh type */
	typedef Mesh mesh_t;

	enum {
		quantity_dimension = Dimension::value
	};

	/** \brief the mesh eleme type vector type */
	typedef typename mesh_t::elem_type_vector_t etv_t;

	/** \brief the field type vector */
	typedef typename tmp::transform<
		etv_t,
		tmp::inserter <
			typename tmp::empty<etv_t>::type,
			tmp::push_back<tmp::_1, tmp::_2>
		>,
		field_view<tmp::_1, FieldOption, Dimension>
	>::type field_type_vector_t;

	/** \brief the iterator class traversing a field subvector */
	template <class field_t>
	struct iterator : casted_iterator<
		typename mesh_elem_iterator_t<typename field_t::elem_t>::type,
		field_t
	> {};
};


/**
* \brief implementation class of NiHu::function_space_view
* \tparam Mesh the underlying Mesh type
* \tparam FieldOption determines how the field is generated from the mesh
*/
template<class Mesh, class FieldOption, class Dimension>
class function_space_impl<function_space_view<Mesh, FieldOption, Dimension> > :
	public Mesh
{
public:
	/** \brief the traits class */
	typedef function_space_traits< function_space_view<Mesh, FieldOption, Dimension> > traits_t;

private:
	/** \brief specialisation of number of dofs for the constant case */
	size_t get_num_dofs_impl(field_option::constant) const
	{
		return Mesh::get_num_elements() * traits_t::quantity_dimension;
	}

	/** \brief specialisation of number of dofs for the isoparametric case */
	size_t get_num_dofs_impl(field_option::isoparametric) const
	{
		return Mesh::get_num_points() * traits_t::quantity_dimension;
	}

public:
	/** \brief return mesh reference */
	Mesh const &get_mesh() const
	{
		return static_cast<Mesh const &> (*this);
	}

	/**
	* \brief first field of given element type
	* \tparam ElemType the element type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_begin() const
	{
		return Mesh::template begin<typename FieldType::elem_t>();
	}

	/**
	* \brief last field of given element type
	* \tparam ElemType the element type to access
	*/
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_end() const
	{
		return Mesh::template end<typename FieldType::elem_t>();
	}

	/**
	* \brief return number of degrees of freedom
	* \return number of degrees of freedom
	*/
	size_t get_num_dofs() const
	{
		return get_num_dofs_impl(FieldOption());
	}
};


/**
* \brief A mesh extended with a Field generating option
* \tparam Mesh the underlying Mesh type
* \tparam FieldOption determines how the field is generated from the mesh
*/
template<class Mesh, class FieldOption, class Dimension = _1d>
class function_space_view :
	public function_space_base<function_space_view<Mesh, FieldOption, Dimension> >,
	public function_space_impl<function_space_view<Mesh, FieldOption, Dimension> >
{
public:
	/** \brief the implementation class */
	typedef function_space_impl<function_space_view<Mesh, FieldOption, Dimension> > impl_t;

	using impl_t::get_mesh;
	using impl_t::field_begin;
	using impl_t::field_end;
	using impl_t::get_num_dofs;
};


/** \brief factory function to transform a mesh into a function space
 * \tparam Mesh the mesh type
 * \tparam Option the field generation option
 * \param [in] msh mesh reference
 * \param [in] dim the function space's dimension
 * \return function space view of the mesh
 */
template <class Mesh, class Option, class Dimension = _1d>
function_space_view<Mesh, Option, Dimension> const &
	create_function_space_view(Mesh const &msh, Option, Dimension dim = Dimension())
{
	return static_cast<function_space_view<Mesh, Option, Dimension> const &>(msh);
}

/** \brief factory function to transform a mesh into an isoparametric function space
 * \tparam Mesh the mesh type
 * \param [in] msh mesh reference
 * \param [in] dim the function space's dimension
 * \return isoparametric function space view of the mesh
 */
template <class Mesh, class Dimension = _1d>
function_space_view<Mesh, field_option::isoparametric, Dimension> const &
	isoparametric_view(Mesh const &msh, Dimension dim = Dimension())
{
	return create_function_space_view(msh, field_option::isoparametric(), dim);
}

/** \brief factory function to transform a mesh into a constant function space
 * \tparam Mesh the mesh type
 * \param [in] msh mesh reference
 * \param [in] dim the function space's dimension
 * \return constant function space view of the mesh
 */
template <class Mesh, class Dimension = _1d>
function_space_view<Mesh, field_option::constant, Dimension> const &
	constant_view(Mesh const &msh, Dimension dim = Dimension())
{
	return create_function_space_view(msh, field_option::constant(), dim);
}


// forward declaration
template <class FuncSpace>
class dirac_space;

/** \brief traits class of a NiHu::dirac_space */
template <class FuncSpace>
struct function_space_traits<dirac_space<FuncSpace> >
{
	/** \brief the parent traits */
	typedef function_space_traits<FuncSpace> base_traits;

	/** \brief the mesh type */
	typedef typename base_traits::mesh_t mesh_t;

	/** \brief the field type vector */
	typedef typename tmp::transform<
		typename base_traits::field_type_vector_t,
		tmp::inserter<
			typename tmp::empty<
				typename base_traits::field_type_vector_t
			>::type,
			tmp::push_back<tmp::_1, tmp::_2>
		>,
		dirac_field<tmp::_1>
	>::type field_type_vector_t;

	/** \brief the iterator class traversing a field subvector */
	template <class dirac_field_t>
	struct iterator : casted_iterator<
		typename base_traits::template iterator<
			typename dirac_field_t::field_t
		>::type,
		dirac_field_t,
		field_impl<typename dirac_field_t::field_t>
	> {};
};


/** \brief Dirac-like extension of a function space
* \tparam FuncSpace the function space type to extend into dirac_space
*/
template <class FuncSpace>
class dirac_space :
	public function_space_base<dirac_space<FuncSpace> >,
	public function_space_impl<FuncSpace>
{
public:
	/** \brief the traits class */
	typedef function_space_traits<dirac_space<FuncSpace> > traits_t;
	/** \brief the implementation class */
	typedef function_space_impl<FuncSpace> impl_t;
	
	// explicitly define mesh_t because it would be ambiguous otherwise
	typedef typename traits_t::mesh_t mesh_t;

	using impl_t::get_mesh;
	using impl_t::get_num_dofs;

	/** \brief return begin iterator of a subvector of dirac fields */
	template <class dirac_field_t>
	typename traits_t::template iterator<dirac_field_t>::type
		field_begin() const
	{
		return impl_t::template field_begin<typename dirac_field_t::field_t>();
	}

	/** \brief return end iterator of a subvector of dirac fields */
	template <class dirac_field_t>
	typename traits_t::template iterator<dirac_field_t>::type
		field_end() const
	{
		return impl_t::template field_end<typename dirac_field_t::field_t>();
	}
};


/** \brief factory function to convert a function space into a dirac space
* \tparam FuncSpace the function space to convert
* \param [in] space the function space reference to convert
* \return the converted dirac space reference
*/
template <class FuncSpace>
dirac_space<FuncSpace> const &
	dirac(function_space_base<FuncSpace> const &space)
{
	return static_cast<dirac_space<FuncSpace> const &>(
		static_cast<function_space_impl<FuncSpace> const &>(space.derived()));
}



// forward declaration
template <class FieldTypeVector>
class function_space;

/** \brief metafunction to return the element type vector of a field type vector */
template <class FieldTypeVector>
struct field_2_elem_type_vector
{
	/** \brief helper metafunction to extract the element type of a field */
	template <class field_t>
	struct elemize
	{
		typedef typename field_t::elem_t type;
	};

	/** \brief the element type vector */
	typedef typename tmp::unique<
		typename tmp::transform<
			FieldTypeVector,
			tmp::inserter<
				typename tmp::empty<FieldTypeVector>::type,
				tmp::push_back<tmp::_1, tmp::_2>
			>,
			elemize<tmp::_1>
		>::type
	>::type type;
};


/** \brief Traits class of a function space */
template <class FieldTypeVector>
struct function_space_traits<function_space<FieldTypeVector> >
{
	/** \brief the underlying mesh type */
	typedef mesh<typename field_2_elem_type_vector<FieldTypeVector>::type> mesh_t;

	/** \brief the field type vector */
	typedef FieldTypeVector field_type_vector_t;

	/** \brief the iterator type of a field type subvector */
	template <class field_t>
	struct iterator
	{
		typedef typename eigen_std_vector<field_t>::type::const_iterator type;
	};
};


/** \brief class describing a function space
* \tparam FieldTypeVector compile time vector of fields building the function space
*/
template <class FieldTypeVector>
class function_space_impl<function_space<FieldTypeVector> > :
	public function_space_traits<function_space<FieldTypeVector> >::mesh_t
{
public:
	/** \brief the traits class */
	typedef function_space_traits<function_space<FieldTypeVector> > traits_t;

	/** \brief the field type vector */
	typedef typename traits_t::field_type_vector_t field_type_vector_t;

	/** \brief the underlying mesh type */
	typedef typename traits_t::mesh_t mesh_t;

	/** \brief combine field_type_vector into a BIG heterogeneous std::vector container */
	typedef typename tmp::inherit<
		typename tmp::transform<
			field_type_vector_t,
			tmp::inserter<
				typename tmp::empty<field_type_vector_t>::type,
				tmp::push_back<tmp::_1,tmp::_2>
			>,
			eigen_std_vector<tmp::_1>
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

	/** \brief return underlying mesh reference
	 * \return the mesh
	 */
	mesh_t const &get_mesh() const
	{
		return static_cast<mesh_t const &>(*this);
	}

	/** \brief return begin iterator of a subspace */
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_begin() const
	{
		return m_fields.eigen_std_vector<FieldType>::type::begin();
	}

	/** \brief return end iterator of a subspace */
	template <class FieldType>
	typename traits_t::template iterator<FieldType>::type
		field_end() const
	{
		return m_fields.eigen_std_vector<FieldType>::type::end();
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
		m_fields.eigen_std_vector<field_t>::type::push_back(f.derived());
		return *(m_fields.eigen_std_vector<field_t>::type::rbegin());
	}

	/**
	 * \brief return number of fields
	 * \return number of fields in the function space
	 */
	unsigned get_num_fields() const
	{
		return this->get_num_elements();
	}

	/** \brief return number of dofs */
	unsigned get_num_dofs() const
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
	/** \brief the implementation type */
	typedef function_space_impl<function_space<FieldTypeVector> > impl_t;

	/** \brief mesh type redefined to avoid ambigous lookup */
	typedef typename impl_t::mesh_t mesh_t;

	using impl_t::field_begin;
	using impl_t::field_end;
	using impl_t::get_num_dofs;
	using impl_t::get_mesh;

	/** \brief constructor */
	function_space() :
		impl_t()
	{
	}

	/** \brief constructor from node and field definition matrices */
	template <class node_matrix_t, class field_matrix_t>
	function_space(node_matrix_t const &nodes, field_matrix_t const &fields)
		: impl_t(nodes, fields)
	{
	}
};

/** \brief factory function to create a function space from fields
 * \tparam nodes_t the nodes matrix type
 * \tparam elements_t the elements matrix type
 * \tparam fields_t the field matrix type
 * \param [in] nodes the nodes matrix
 * \param [in] elements the field description matrix
 * \param [in] fields the field tag instances
 */
template <class nodes_t, class elements_t, class...fields_t>
function_space<tmp::vector<typename tag2field<fields_t>::type...> >
	create_function_space(nodes_t const &nodes, elements_t const &elements, fields_t const &...fields)
{
	return function_space<tmp::vector<typename tag2field<fields_t>::type...> >(nodes, elements);
}

/** \brief metafunction determining if argument is function space expression
 * \tparam FuncSpace the function space to test
 */
template <class FuncSpace>
struct is_function_space : std::is_base_of<
	function_space_base<typename std::decay<FuncSpace>::type>,
	typename std::decay<FuncSpace>::type
> {};

}

#endif // FUNCTION_SPACE_HPP_INCLUDED

