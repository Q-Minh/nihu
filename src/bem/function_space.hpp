#ifndef FUNCTION_SPACE_HPP_INCLUDED
#define FUNCTION_SPACE_HPP_INCLUDED

#include "field.hpp"
#include "mesh.hpp"

struct constant;
struct isoparametric;

template <class ElemType, class Method>
struct generate_field;

template <class ElemType>
struct generate_field<ElemType, constant>
{
	typedef ConstantField<ElemType> type;
};

template <class ElemType>
struct generate_field<ElemType, isoparametric>
{
	typedef IsoParametricField<ElemType> type;
};

template<class ElemVector, class FieldGenerationMethod>
class FunctionSpace
{
public:
	typedef ElemVector elem_type_vector_t;
	typedef FieldGenerationMethod generation_method;

	typedef typename transform<
		elem_type_vector_t,
		inserter<tiny<>, push_back<_1, _2> >,
		generate_field<_1, generation_method>
	>::type field_type_vector_t;

	typedef typename inherit<
		typename transform<
		field_type_vector_t,
		inserter<tiny<>, push_back<_1,_2> >,
		vectorize<_1>
		>::type
	>::type field_container_t;

	FunctionSpace(Mesh<elem_type_vector_t> const &mesh)
	{
	}

protected:
	field_container_t fields;
};

#endif
