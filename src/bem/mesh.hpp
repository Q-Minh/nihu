/**
 * \file mesh.hpp
 * \author Peter fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 * \brief Decalaration of class Mesh
 */
#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include "../tmp/integer.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/algorithm.hpp"
#include "../tmp/control.hpp"

#include "element.hpp"

#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >

#include <iostream>

template <class ElemTypeVector>
class Mesh
{
public:
	/** \brief define template parameter as nested type */
	typedef ElemTypeVector elem_type_vector_t;

	/** \brief the type of the mesh class itself */
	typedef Mesh<elem_type_vector_t> mesh_t;

	/** \brief metafunction to convert T into std::vector<T> */
	template <class T>
	struct vectorize { typedef EIGENSTDVECTOR(T) type; };

	/** \brief combine elem_vector into a BIG heterogeneous std::vector container */
	typedef typename inherit<
		typename transform<
		elem_type_vector_t,
		inserter<tiny<>, push_back<_1,_2> >,
		vectorize<_1>
		>::type
	>::type elem_container_t;

	/** \brief type of a nodal vector */
	typedef typename deref<
		typename begin<elem_type_vector_t>::type
	>::type::x_t x_t;

	/** \brief type of a nodal coordinate */
	typedef typename x_t::Scalar scalar_t;

	/** \brief number of dimensions */
	static unsigned const nDim = x_t::SizeAtCompileTime;

protected:
	EIGENSTDVECTOR(x_t) nodes;	/**< \brief nodal coordinates */
	elem_container_t elements;	/**< \brief element geometries (BIG heterogeneous container) */

	/** \brief wrapper class of member function print_coords used by call_each */ 
	template <class ElemType>
	struct printCoords
	{
		struct type
		{
			void operator() (mesh_t &m)
			{
				m.print_coords<ElemType>();
			}
		};
	};

	template <class elem_t>
	struct elem_adder
	{
		struct type
		{
			bool operator() (unsigned const input[], mesh_t &m)
			{
				if (input[0] == elem_t::domain_t::id)
				{
					// construct element
					typename elem_t::nodes_t nodes;
					typename elem_t::coords_t coords;
					for (unsigned i = 0; i < elem_t::num_nodes; ++i)
					{
						nodes[i] = input[i+1];
						coords.row(i) = m.nodes[nodes[i]];
					}
					m.push_element(elem_t(0, nodes, coords));
				}
				return false;
			}
		};
	};

public:
	/**
	 * \brief return begin iterator of the elements
	 */
	template <class ElemType>
	typename EIGENSTDVECTOR(ElemType)::const_iterator elembegin(void) const
	{
		return elements.EIGENSTDVECTOR(ElemType)::begin();
	}
	
	/**
	 * \brief return end iterator of the elements
	 */
	template <class ElemType>
	typename EIGENSTDVECTOR(ElemType)::const_iterator elemend(void) const
	{
		return elements.EIGENSTDVECTOR(ElemType)::end();
	}
	
	/**
	 * \biref print element coordinates
	 * \tparam ElemType the eleme type that is processed
	 */
	template <class ElemType>
	void print_coords()
	{
		std::cout << ElemType::domain_t::id << " " << ElemType::num_nodes << std::endl;
		std::for_each(
			elements.EIGENSTDVECTOR(ElemType)::begin(),
			elements.EIGENSTDVECTOR(ElemType)::end(),
			[] (ElemType const &e) { std::cout << e.get_coords() << std::endl << std::endl; }
		);
	}

	/**
	 * \brief add a new element to the mesh
	 * \param input array of unsigned values.
	 * input[0] is the elem ID, the subsequent elements are the nodal indices in the mesh
	 * \returns true if the element is inserted into the mesh
	 */
	bool add_elem(unsigned const input[])
	{
		return tmp::call_until<
			elem_type_vector_t,
			elem_adder<_1>,
			unsigned const*,
			mesh_t &
		>(input, *this);
	}

	void add_node(scalar_t input[])
	{
		x_t c;
		typedef typename x_t::Index index_t;
		for (index_t i = 0; i < nDim; ++i)
			c[i] = input[i];
		nodes.push_back(c);
	}

	template <class elem_t>
	void push_element(elem_t const &e)
	{
		elements.EIGENSTDVECTOR(elem_t)::push_back(e);
	}

	void print_elements(void)
	{
		tmp::call_each<
			elem_type_vector_t,
			printCoords<_1>,
			mesh_t &
		>(*this);
	}
};

#endif

