#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>

template <unsigned nDim>
class Coord : public std::array<double, nDim>
{
	typedef Coord type; /* self-returning structure */
};

class tria_tag;		/* triangular element */
class lin_quad_tag;	/* quad element with constant Jacobian (parallelogram) */
class plane_quad_tag;	/* quad element with constant normal */
class quad_tag;	/* general quad element */

/** \brief metafunction returning number of element nodes */
template <class tag> class nNodes;
template<> struct nNodes<tria_tag> : int_<3> {};
template<> struct nNodes<lin_quad_tag> : int_<4> {};
template<> struct nNodes<plane_quad_tag> : int_<4> {};
template<> struct nNodes<quad_tag> : int_<4> {};

template <class ElemTag>
class Elem : public std::array<unsigned, nNodes<ElemTag>::value>
{
public:
	typedef ElemTag elem_tag;
	typedef Elem type;	/* self-returning structure */
	static unsigned const nNod = nNodes<elem_tag>::value;

	template <class nodeIterator>
	bool build(unsigned input[], nodeIterator iter)
	{
		/*
		typedef nodeIterator::value_type node_type;
		std::array<node_type,nNod> nodes;
		*/

		if (ElemTypeRecogniser<elem_tag>::eval(input, iter))
		{
			for (unsigned i = 0; i < nNod; ++i)
				(*this)[i] = *(input+i+1);
			return true;
		}
		return false;
	}
};

template <class elem_tag>
class ElemTypeRecogniser;

template <>
class ElemTypeRecogniser<tria_tag>
{
public:
	template <class nodeIterator>
	static bool eval(unsigned input[], nodeIterator)
	{
		return (input[0] == nNodes<tria_tag>::value);
	}
};

template <>
class ElemTypeRecogniser<lin_quad_tag>
{
public:
	template <class nodeIterator>
	static bool eval(unsigned input[], nodeIterator)
	{
		return (input[0] == nNodes<lin_quad_tag>::value);
	}
};

template <>
class ElemTypeRecogniser<quad_tag>
{
public:
	template <class nodeIterator>
	static bool eval(unsigned input[], nodeIterator)
	{
		return (input[0] == nNodes<quad_tag>::value);
	}
};

typedef Elem<tria_tag> TriaElem;
typedef Elem<lin_quad_tag> LinQuadElem;
typedef Elem<plane_quad_tag> PlaneQuadElem;
typedef Elem<quad_tag> QuadElem;

#endif
