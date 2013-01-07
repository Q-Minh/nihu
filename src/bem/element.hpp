#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <array>

template <unsigned nDim>
class Coord : public std::array<double, nDim>
{
	typedef Coord type; /* self-returning structure */
};


class tria_tag;
class quad_tag;

/** \brief metafunction returning number of element nodes */
template <class tag> class nNodes;
template<> struct nNodes<quad_tag> : int_<3> {};
template<> struct nNodes<tria_tag> : int_<4> {};

template <class Elem_tag>
class Elem : public std::array<unsigned, nNodes<Elem_tag>::value>
{
public:
	typedef Elem_tag elem_tag;
	typedef Elem type;	/* self-returning structure */
	static unsigned const nNod = nNodes<Elem_tag>::value;

	bool build(unsigned input[])
	{
		if (input[0] == nNod)
		{
			for (unsigned i = 0; i < nNod; ++i)
				(*this)[i] = *(input+i+1);
			return true;
		}
		return false;
	}
};

typedef Elem<quad_tag> QuadElem;
typedef Elem<tria_tag> TriaElem;

#endif

