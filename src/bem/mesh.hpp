#ifndef MESH_HPP
#define MESH_HPP

#include "../tmp/integer.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/algorithm.hpp"
#include "../tmp/call_each.hpp"

#include "element.hpp"

#include <vector>

/** \brief metafunction to convert T into vector<T> */
template <class T>
struct vectorize { typedef std::vector<T> type; };


template <unsigned nDim, class ElemTypes>
class Mesh
{
private:
	/** \brief transform sequence of ElemType into sequence of vector<ElemType> T */
	typedef typename transform<
		typename begin<ElemTypes>::type,
		typename end<ElemTypes>::type,
		inserter<tiny<>, push_back<_1,_2> >,
		vectorize<_1>
	>::type elem_vectors;

	/** \brief combine vectors into a big heterogeneous vector container T */
	typedef typename inherit<
		typename begin<elem_vectors>::type,
		typename end<elem_vectors>::type
	>::type ElemVector;

	std::vector<Coord<3> > nodes;	/**< \brief nodal coordinates */
	ElemVector elements;			/**< element nodes (heterogeneous container) */

	/** \brief functor template used by call_each to build add_elem member function's for loop */
	template <class ElemType>
	struct elem_adder
	{
		void operator() (unsigned input[], ElemVector &v)
		{
			ElemType e;
			if (e.build(input))
				v.std::vector<ElemType>::push_back(e);
		}
	};

public:

	void add_elem(unsigned input[])
	{
		/* build class that subsequently invokes elem_adder<T>() for each element type T */
		typedef typename call_each<
			typename begin<ElemTypes>::type,
			typename end<ElemTypes>::type,
			elem_adder<_1>,
			unsigned*,
			ElemVector &
		>::type b;

		/* here are the calls invoked */
		b::apply(input, elements);
	}
};

#endif

