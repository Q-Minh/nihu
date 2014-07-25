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
 * \file mesh.hpp
 * \ingroup funcspace
 * \brief Declaration of class Mesh
 */
#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include "../tmp/integer.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/algorithm.hpp"
#include "../tmp/control.hpp"
#include "tmp/vector.hpp"

#include "element.hpp"
#include "../util/eigen_utils.hpp"

/** \brief metafunction returning the iterator that traverses the homogeneous element vector of specified element type
 * \tparam ElemType the specified element type
 */
template <class ElemType>
struct mesh_elem_iterator_t
{
	/** \brief the iterator that traverses a submesh's element container */
	typedef typename EigenStdVector<ElemType>::type::const_iterator type;
};


/**
 * \brief container class for field points
 * \tparam xType type of a field point
 */
template <class xType>
class field_points
{
public:
	typedef xType x_t; /**< \brief template parameter as nested type */
	static unsigned const nDim = x_t::SizeAtCompileTime; /**< \brief number of dimensions */
	typedef typename EigenStdVector<x_t>::type::const_iterator iterator_t;	/**< \brief node iterator type */

	/**
	 * \brief add a point to the field point mesh
	 * \param p the point to be added
	 */
	void add_point(x_t const &p)
	{
		points.push_back(p);
	}

	/**
	 * \brief return number of points
	 * \return number of points
	 */
	unsigned get_num_points(void) const
	{
		return points.size();
	}

protected:
	typename EigenStdVector<x_t>::type points;	/**< \brief nodal coordinates */
};

/** \brief metafunction computing the first element's x_t in a vector of elements */
template <class ElemTypeVector>
struct first_elements_x_type
{
	typedef typename tmp::deref<
		typename tmp::begin<ElemTypeVector>::type
	>::type::x_t type;
};

/**
 * \brief container class for a mesh
 * \tparam ElemTypeVector compile time vector of the contained element types
 */
template <class ElemTypeVector>
class mesh :
	public field_points<typename first_elements_x_type<ElemTypeVector>::type>
{
public:
	/** \brief define template parameter as nested type */
	typedef ElemTypeVector elem_type_vector_t;

	/** \brief type of base class */
	typedef field_points<typename first_elements_x_type<ElemTypeVector>::type> base_t;

	/** \brief number of dimensions of the mesh */
	static unsigned const nDim = base_t::nDim;

	/** \brief combine elem_vector into a BIG heterogeneous std::vector container */
	typedef typename tmp::inherit<
		typename tmp::transform<
			elem_type_vector_t,
			tmp::inserter<
				typename tmp::empty<elem_type_vector_t>::type,
				tmp::push_back<tmp::_1,tmp::_2>
			>,
			EigenStdVector<tmp::_1>
		>::type
	>::type elem_container_t;

	/** \brief type of a nodal vector */
	typedef typename base_t::x_t x_t;

	/** \brief type of a nodal coordinate */
	typedef typename x_t::Scalar scalar_t;

	template <class ElemType>
	struct elem_iterator_t : mesh_elem_iterator_t<ElemType> {};


protected:
	template <class elem_t>
	struct elem_adder { struct type	{
		/**
		 * \brief add an element of given type to the mesh
		 * \param m the mesh to extend
		 * \param input array containing element node indices
		 */
		bool operator() (unsigned const input[], mesh &m)
		{
			if (input[0] == elem_t::id)
			{
				// construct element
				typename elem_t::nodes_t nodes;
				typename elem_t::coords_t coords;
				for (unsigned i = 0; i < elem_t::num_nodes; ++i)
				{
					nodes[i] = input[i+1];
					coords.col(i) = m.points[nodes[i]];
				}
				m.push_element(elem_t(coords, m.m_num_elements++, nodes));
                return true;
			}
			return false;
		}
	};};

public:
	mesh() : m_num_elements(0)
	{
	}

	/**
	 * \brief build the mesh from mesh description matrices
	 * \tparam node_t type of nodes matrix
	 * \tparam elem_t type of elements matrix
	 * \param [in] nodes matrix of nodal coordinates
	 * \param [in] elements matrix of element node indices
	 */
	template <class node_t, class elem_t>
	mesh(node_t const &nodes, elem_t const &elements) : m_num_elements(0)
	{
		unsigned const N_MAX_ELEM = 1+9;

		for (unsigned i = 0; i < unsigned(nodes.rows()); ++i)
		{
			double c[nDim];
			for (unsigned j = 0; j < nDim; ++j)
				c[j] = nodes(i,j);
			add_node(c);
		}
		for (unsigned i = 0; i < unsigned(elements.rows()); ++i)
		{
			unsigned e[N_MAX_ELEM] = {0};
			for (unsigned j = 0; j < unsigned(elements.cols()); ++j)
				e[j] = unsigned(elements(i,j));
			add_elem(e);
		}
	}

	/**
	 * \brief return begin iterator of the elements
	 */
	template <class ElemType>
	typename elem_iterator_t<ElemType>::type begin(void) const
	{
		return m_elements.EigenStdVector<ElemType>::type::begin();
	}

	/**
	 * \brief return end iterator of the elements
	 */
	template <class ElemType>
	typename elem_iterator_t<ElemType>::type end(void) const
	{
		return m_elements.EigenStdVector<ElemType>::type::end();
	}

	/**
	 * \brief add a new element to the mesh
	 * \param input array of unsigned values.
	 * input[0] is the elem ID, the subsequent elements are the nodal indices in the mesh
	 * \return true if the element is inserted into the mesh
	 */
	bool add_elem(unsigned const input[])
	{
		return tmp::call_until<
			elem_type_vector_t,
			elem_adder<tmp::_1>,
			unsigned const*,
			mesh &
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
	template <class elem_t>
	elem_t const &push_element(element_base<elem_t> const &e)
	{
		m_elements.EigenStdVector<elem_t>::type::push_back(e.derived());
		return *(m_elements.EigenStdVector<elem_t>::type::rbegin());
	}

	/**
	 * \brief return number of elements
	 * \return number of elements in the mesh
	 */
	unsigned get_num_elements(void) const
	{
		return m_num_elements;
	}

protected:
	elem_container_t m_elements;	/**< \brief element geometries (BIG heterogeneous container) */
	unsigned m_num_elements;	/**< \brief total number of elements in the mesh */
};


/** \brief factory function to create a mesh from nodes and elements matrices
 * \tparam nodes_t type of the node definition matrix
 * \tparam elements_t type of the element definition matrix
 * \tparam Args element types of the mesh
 * \param [in] nodes the node definition matrix instance
 * \param [in] elements the element definition matrix instance
 * \return a mesh consisting of given element types
 */
template <class nodes_t, class elements_t, class...Args>
mesh<tmp::vector<typename tag2element<Args>::type...> >
	create_mesh(nodes_t const &nodes, elements_t const &elements, Args...)
{
	return mesh<tmp::vector<typename tag2element<Args>::type...> >(nodes, elements);
}


template <class Mesh, class Elem>
class homogeneous_submesh
{
public:
	typedef Mesh mesh_t;
	typedef Elem elem_t;

	homogeneous_submesh(mesh_t const &parent) : m_parent(parent)
	{
	}

	/** \brief return begin iterator of the elements */
	typename mesh_t::template elem_iterator_t<elem_t>::type begin(void) const
	{
		return m_parent.template begin<elem_t>();
	}

	/** \brief return end iterator of the elements */
	typename mesh_t::template elem_iterator_t<elem_t>::type end(void) const
	{
		return m_parent.template end<elem_t>();
	}

private:
	mesh_t const &m_parent;
};

/** \brief factory function to create a homogeneous submesh from a mesh
 * \tparam Mesh the mesh type
 * \tparam Tag the element tag type
 * \param [in] mesh the parent mesh
 */
template <class Mesh, class Tag>
homogeneous_submesh<Mesh, typename tag2element<Tag>::type>
create_homogeneous_submesh(Mesh const &mesh, Tag)
{
	return homogeneous_submesh<Mesh, typename tag2element<Tag>::type >(mesh);
}

#endif

