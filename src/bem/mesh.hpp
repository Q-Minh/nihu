/**
 * \file mesh.hpp
 * \author Peter fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 * \brief Declaration of class Mesh
 */
#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include "../tmp/integer.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/algorithm.hpp"
#include "../tmp/control.hpp"

#include "element.hpp"

#include <Eigen/StdVector>
/**
 * \brief macro declaring an Eigen std::vector type with the appropriate allocator
 */
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
	typedef xType x_t; /**< \brief template parameter as nested type */
	static unsigned const nDim = x_t::SizeAtCompileTime; /**< \brief number of dimensions */
	typedef typename EIGENSTDVECTOR(x_t)::const_iterator iterator_t;	/**< \brief node iterator type */

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

/**
 * \brief container class for a mesh
 * \tparam ElemTypeVector compile time vector of the contained element types
 */
template <class ElemTypeVector>
class Mesh : public field_points<typename first_elements_x_type<ElemTypeVector>::type>
{
public:
	/** \brief define template parameter as nested type */
	typedef ElemTypeVector elem_type_vector_t;

	/** \brief type of base class */
	typedef field_points<typename first_elements_x_type<ElemTypeVector>::type> base;

	static unsigned const nDim = base::nDim;	/**< \brief number of dimensions of the mesh */

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

	template <class ElemType>
	struct elem_iterator_t
	{
		typedef typename EIGENSTDVECTOR(ElemType)::const_iterator type;
	};

protected:
	elem_container_t elements;	/**< \brief element geometries (BIG heterogeneous container) */
	unsigned num_elements;	/**< \brief total number of elements in the mesh */

	template <class elem_t>
	struct elem_adder { struct type	{
		/**
		 * \brief add an element of given type to the mesh
		 * \param m the mesh to extend
		 * \param input array containing element node indices
		 */
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
	};};

public:
	Mesh() : num_elements(0)
	{
	}

	/**
	 * \brief build the mesh from MATLAB matrices
	 * \tparam N_MAX_ELEM the maximal width of the elements matrix
	 * \tparam type of nodes matrix
	 * \tparam type of elements matrix
	 * \param [in] nodes matrix of nodal coordinates
	 * \param [in] elements matrix of element node indices
	 */
	template <unsigned N_MAX_ELEM, class node_t, class elem_t>
	void build_from_mex(node_t const &nodes, elem_t const &elements)
	{
		double c[nDim];

		for (size_t i = 0; i < nodes.rows(); ++i)
		{
			for (unsigned j = 0; j < nDim; ++j)
				c[j] = nodes(i,j);
			add_node(c);
		}
		unsigned e[N_MAX_ELEM];
		for (size_t i = 0; i < elements.rows(); ++i)
		{
			for (unsigned j = 0; j < N_MAX_ELEM; ++j)
				e[j] = (unsigned)elements(i,j);
			add_elem(e);
		}
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
	void push_element(elem_t const &e)
	{
		elements.EIGENSTDVECTOR(elem_t)::push_back(e);
		++num_elements;
	}

	/**
	 * \brief return number of elements
	 * \return number of elements in the mesh
	 */
	unsigned get_num_elements(void) const
	{
		return num_elements;
	}
};

#endif

