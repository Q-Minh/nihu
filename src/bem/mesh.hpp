#ifndef MESH_HPP
#define MESH_HPP

#include "../tmp/integer.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/algorithm.hpp"
#include "../tmp/call_each.hpp"

#include "element.hpp"
#include "node.hpp"

#include <vector>

/** \brief metafunction to convert T into vector<T> */
template <class T>
struct vectorize { typedef std::vector<T> type; };


template <int nDim, class ElemTypes>
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

	std::vector<Coord<nDim> > nodes;	/**< \brief nodal coordinates */

	/** \brief functor template used by call_each to build add_elem member function's for loop */
	template <class ElemType>
	struct elem_adder
	{
		struct type {
			bool operator() (int input[], Mesh<nDim, ElemTypes> &m)
			{
				ElemType e;
				if (e.build(input, m.nodes.begin()))
				{
					m.elements.std::vector<ElemType>::push_back(e);
					return true;
				}

				return false;
			}
		};
	};

public:

	ElemVector elements;			/**< element nodes (heterogeneous container) */
	
	bool add_elem(int input[])
	{
		return call_until<ElemTypes, elem_adder<_1> >(input, *this);
	}

	void add_node(double input[])
	{
		Coord<nDim> c;
		c << input[0], input[1], input[2];
		nodes.push_back(c);
	}
};

#endif


