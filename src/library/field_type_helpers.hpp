#ifndef NIHU_FIELD_TYPE_HELPERS_HPP_INCLUDED
#define NIHU_FIELD_TYPE_HELPERS_HPP_INCLUDED

#include "lib_shape.hpp"

namespace NiHu
{

template <class Field>
struct is_constant_line : std::integral_constant<bool,
	std::is_same<typename Field::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename Field::nset_t, line_0_shape_set>::value
> {};

template <class Field>
struct is_linear_line : std::integral_constant<bool,
	std::is_same<typename Field::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename Field::nset_t, line_1_shape_set>::value
> {};

template <class Field>
struct is_constant_tria : std::integral_constant<bool,
	std::is_same<typename Field::elem_t::lset_t, tria_1_shape_set>::value &&
	std::is_same<typename Field::nset_t, tria_0_shape_set>::value
> {};

} // end of namespace NiHu


#endif