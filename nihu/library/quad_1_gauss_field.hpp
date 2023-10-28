#ifndef QUAD_1_GAUSS_FIELD_HPP_INCLUDED
#define QUAD_1_GAUSS_FIELD_HPP_INCLUDED

#include "quad_1_gauss_shape_set.hpp"
#include "../core/field.hpp"
#include "../util/type2tag.hpp"
#include "lib_element.hpp"

// define traits of the new field type
namespace NiHu
{

// define the new field type
typedef field<
	quad_1_elem,			// over a linear quad elem
	quad_1_gauss_shape_set,	// using the new shape set
	field_dimension::_1d 	// scalar field
> quad_1_gauss_field;

namespace field_traits
{
	template <>
	struct id<quad_1_gauss_field> { enum {value = 666}; };
} // end of namespace field_traits

// define a tag to the new type
typedef type2tag<quad_1_gauss_field> quad_1_gauss_field_tag;

} // end of namespace NiHu

#endif // QUAD_1_GAUSS_FIELD_HPP_INCLUDED
