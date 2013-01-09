#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <Eigen/Dense>
using Eigen::Matrix;

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

template <class elem_tag>
class ElemTypeRecogniser;


template <class ElemTag>
class Elem : public Matrix<int, nNodes<ElemTag>::value, 1>
{
public:
	typedef ElemTag elem_tag;
	typedef Elem type;	/* self-returning structure */
	static int const nNod = nNodes<elem_tag>::value;

	template <class nodeIterator>
	bool build(int input[], nodeIterator iter)
	{
		if (ElemTypeRecogniser<elem_tag>::eval(input, iter))
		{
			for (int i = 0; i < nNod; ++i)
				(*this)[i] = *(input+i+1);
			return true;
		}
		return false;
	}
};

template <>
class ElemTypeRecogniser<tria_tag>
{
public:
	template <class nodeIterator>
	static bool eval(int input[], nodeIterator)
	{
		return (input[0] == nNodes<tria_tag>::value);
	}
};

template <>
class ElemTypeRecogniser<lin_quad_tag>
{
public:
	template <class nodeIterator>
	static bool eval(int input[], nodeIterator iter)
	{
		typedef typename nodeIterator::value_type node_type;
		
		if (input[0] == nNodes<lin_quad_tag>::value)
		{
			node_type c0 = *(iter+input[1]);
			node_type c1 = *(iter+input[2]);
			node_type c2 = *(iter+input[3]);
			node_type c3 = *(iter+input[4]);
			
			return ((c1-c0) + (c3-c0) - (c2-c0)).norm() < 1e-3;
		}
		else
			return false;
	}
};

template <>
class ElemTypeRecogniser<plane_quad_tag>
{
public:
	template <class nodeIterator>
	static bool eval(int input[], nodeIterator iter)
	{
		typedef typename nodeIterator::value_type node_type;
		
		if (input[0] == nNodes<plane_quad_tag>::value)
		{
			node_type c0 = *(iter+input[1]);
			node_type c1 = *(iter+input[2]);
			node_type c2 = *(iter+input[3]);
			node_type c3 = *(iter+input[4]);
			
			return ((c1-c0).cross(c2-c0)).dot(c3-c0) < 1e-3;
		}
		else
			return false;
	}
};

template <>
class ElemTypeRecogniser<quad_tag>
{
public:
	template <class nodeIterator>
	static bool eval(int input[], nodeIterator)
	{
		return input[0] == nNodes<lin_quad_tag>::value;
	}
};

typedef Elem<tria_tag> TriaElem;
typedef Elem<lin_quad_tag> LinQuadElem;
typedef Elem<plane_quad_tag> PlaneQuadElem;
typedef Elem<quad_tag> QuadElem;

#endif

