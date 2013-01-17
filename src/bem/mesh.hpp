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
	ElemVector elements;		/**< element nodes (heterogeneous container) */

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
				else
					return false;
			}
		};
	};

public:

	bool add_elem(unsigned input[])
	{
		switch (input[0])
		{
		case 3:
			{
				tria_1_elem::coords_t coords;
				for (unsigned i = 0; i < tria_1_elem::num_nodes; ++i)
					coords.row(i) = nodes[input[1+i]];
				elements.std::vector<tria_1_elem, Eigen::aligned_allocator<tria_1_elem> >::push_back(tria_1_elem(coords));
			}
			break;
		case 4:
			{
				quad_1_elem::coords_t coords;
				for (unsigned i = 0; i < quad_1_elem::num_nodes; ++i)
					coords.row(i) = nodes[input[1+i]];

				typedef tiny<parallelogram_elem, quad_1_elem> smallElemTypes;

				call_until<
					smallElemTypes,
					elem_adder<quad_1_elem, _1>,
					quad_1_elem::coords_t const &, Mesh<ElemTypes> &
				>(coords, *this);
			}
			break;
		default:
			return false;
		}
		return true;
	}

	void add_node(scalar_t input[])
	{
		x_t c;
		typedef x_t::Index index_t;
		for (index_t i = 0; i < nDim; ++i)
			c[i] = input[i];
		nodes.push_back(c);
	}
};

#endif

