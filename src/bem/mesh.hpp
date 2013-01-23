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

/**
 * \brief container class for field points
 * \tparam xType type of a field point
 */
template <class xType>
class field_points
{
public:
	typedef xType x_t;
	typedef typename EIGENSTDVECTOR(x_t)::const_iterator iterator_t;

	iterator_t begin(void) const
	{
		return points.begin();
	}

	iterator_t end(void) const
	{
		return points.end();
	}

	void add_point(x_t const &p)
	{
		points.push_back(p);
	}

	unsigned get_num_points(void) const
	{
		return points.size();
	}

protected:
	EIGENSTDVECTOR(x_t) points;	/**< \brief nodal coordinates */
};

/** \brief metafunction computing the first element's x_t in a vector of elements */
template <class ElemTypeVector>
struct first_elements_x_type
{
	typedef typename tmp::deref<
		typename tmp::begin<ElemTypeVector>::type
	>::type::x_t type;
};

template <class ElemTypeVector>
class Mesh : public field_points<typename first_elements_x_type<ElemTypeVector>::type>
{
public:
	/** \brief define template parameter as nested type */
	typedef ElemTypeVector elem_type_vector_t;

	/** \brief type of base class */
	typedef field_points<typename first_elements_x_type<ElemTypeVector>::type> base;

	/** \brief the type of the mesh class itself */
	typedef Mesh<elem_type_vector_t> mesh_t;

	/** \brief metafunction to convert T into std::vector<T> */
	template <class T>
	struct vectorize { typedef EIGENSTDVECTOR(T) type; };

	/** \brief combine elem_vector into a BIG heterogeneous std::vector container */
	typedef typename tmp::inherit<
		typename tmp::transform<
		elem_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1,tmp::_2> >,
		vectorize<tmp::_1>
		>::type
	>::type elem_container_t;

	/** \brief type of a nodal vector */
	typedef typename base::x_t x_t;

	/** \brief type of a nodal coordinate */
	typedef typename x_t::Scalar scalar_t;

	/** \brief number of dimensions */
	static unsigned const nDim = x_t::SizeAtCompileTime;

	template <class ElemType>
	struct elem_iterator_t
	{
		typedef typename EIGENSTDVECTOR(ElemType)::const_iterator type;
	};

protected:
	elem_container_t elements;	/**< \brief element geometries (BIG heterogeneous container) */
	unsigned num_elements;

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
						coords.row(i) = m.points[nodes[i]];
					}
					m.push_element(elem_t(0, nodes, coords));
				}
				return false;
			}
		};
	};

public:
	Mesh() : num_elements(0)
	{
	}

	/**
	 * \brief return begin iterator of the elements
	 */
	template <class ElemType>
	typename elem_iterator_t<ElemType>::type begin(void) const
	{
		return elements.EIGENSTDVECTOR(ElemType)::begin();
	}

	/**
	 * \brief return end iterator of the elements
	 */
	template <class ElemType>
	typename elem_iterator_t<ElemType>::type end(void) const
	{
		return elements.EIGENSTDVECTOR(ElemType)::end();
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
			elem_adder<tmp::_1>,
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
		add_point(c);
	}

	template <class elem_t>
	void push_element(elem_t const &e)
	{
		elements.EIGENSTDVECTOR(elem_t)::push_back(e);
		++num_elements;
	}

	unsigned get_num_elements(void) const
	{
		return num_elements;
	}
};

#endif

