#ifndef LINE_1_GAUSS_FIELD_HPP_INCLUDED
#define LINE_1_GAUSS_FIELD_HPP_INCLUDED

#include "line_1_gauss_shape_set.hpp"
#include "core/field.hpp"
#include "lib_element.hpp"

// define traits of the new field type
namespace NiHu
{

// define the new field type
typedef field<
	line_1_elem,			// over a linear line elem
	line_1_gauss_shape_set,	// using the new shape set
	_1d 					// scalar field
> line_1_gauss_field;

namespace field_traits
{
	template <>
	struct id<line_1_gauss_field> { enum {value = 222}; };
} // end of namespace field_traits

// define a tag to the new type
struct line_1_gauss_field_tag {};

template <>
struct tag2field<line_1_gauss_field_tag> :
	line_1_gauss_field {};
	
} // end of namespace NiHu

#endif // LINE_1_GAUSS_FIELD_HPP_INCLUDED
