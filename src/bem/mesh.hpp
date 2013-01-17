#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include "../tmp/integer.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/algorithm.hpp"
#include "../tmp/call_each.hpp"

#include "element.hpp"

#include <Eigen/StdVector>

/** \brief metafunction to convert T into vector<T> */
template <class T>
struct vectorize { typedef std::vector<T, Eigen::aligned_allocator<T> > type; };


template <class ElemTypes>
class Mesh
{
private:
	/** \brief combine ElemTypes into a big heterogeneous vector container */
	typedef typename inherit<
		typename transform<
		ElemTypes,
		inserter<tiny<>, push_back<_1,_2> >,
		vectorize<_1>
		>::type
	>::type ElemVector;

	/** \brief type of a nodal vector */
	typedef typename deref<  typename begin<ElemTypes>::type  >::type::x_t x_t;
	/** \brief type of a nodal coordinate */
	typedef typename x_t::Scalar scalar_t;
	/** \brief number of dimensions */
	static unsigned const nDim = x_t::SizeAtCompileTime;
	
	std::vector<x_t> nodes;		/**< \brief nodal coordinates */
	ElemVector elements;		/**< \brief element geometries (heterogeneous container) */
	
	template <class RefElem, class ElemType>
	struct elem_adder
	{
		typedef typename RefElem::coords_t coords_t;
		typedef shape_set_converter<typename RefElem::lset_t, typename ElemType::lset_t, nDim> converter_t;

		struct type {
			bool operator() (coords_t const &coords, Mesh<ElemTypes> &m)
			{
				if (converter_t::eval(coords))
				{
					m.elements.std::vector<ElemType, Eigen::aligned_allocator<ElemType> >::push_back(
						ElemType(converter_t::get_coords())
						);
					return true;
				}
				return false;
			}
		};
	};

	template <class smallElemTypes>
	struct elem_adder_outer
	{
		struct type {
			bool operator() (unsigned const input[], Mesh<ElemTypes> &m)
			{
				typedef typename deref<
					typename prev<
						typename end<smallElemTypes>::type
					>::type
				>::type ref_elem_t;

				if (input[0] == ref_elem_t::domain_t::id)
				{
					typename ref_elem_t::coords_t coords;
					for (unsigned i = 0; i < ref_elem_t::num_nodes; ++i)
						coords.row(i) = m.nodes[input[1+i]];

					return call_until<
						smallElemTypes,
						elem_adder<ref_elem_t, _1>,
						typename ref_elem_t::coords_t const &, Mesh<ElemTypes> &
					>(coords, m);
				}
				return false;
			}
		};
	};


public:

	bool add_elem(unsigned const input[])
	{
		typedef tiny< tiny<tria_1_elem>, tiny<parallelogram_elem, quad_1_elem> > elemArray;
		return call_until<
			elemArray,
			elem_adder_outer<_1>,
			unsigned const*, Mesh<ElemTypes> &
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
};

#endif

