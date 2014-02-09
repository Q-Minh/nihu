#ifndef LIB_ELEMENT_HPP_INCLUDED
#define LIB_ELEMENT_HPP_INCLUDED

#include "../core/element.hpp"
#include "lib_shape.hpp"

typedef general_surface_element<line_1_shape_set, space_2d::scalar_t> line_1_elem;
typedef general_surface_element<line_2_shape_set, space_2d::scalar_t> line_2_elem;
typedef general_surface_element<tria_1_shape_set, space_3d::scalar_t> tria_1_elem;
typedef general_surface_element<tria_2_shape_set, space_3d::scalar_t> tria_2_elem;
typedef general_surface_element<quad_1_shape_set, space_3d::scalar_t> quad_1_elem;
typedef general_surface_element<quad_2_shape_set, space_3d::scalar_t> quad_2_elem;
typedef general_surface_element<quad_28_shape_set, space_3d::scalar_t> quad_28_elem;

struct line_1_tag {};
template <> struct tag2element<line_1_tag> : line_1_elem {};
struct line_2_tag {};
template <> struct tag2element<line_2_tag> : line_2_elem {};


struct tria_1_tag {};
template <> struct tag2element<tria_1_tag> : tria_1_elem {};
struct tria_2_tag {};
template <> struct tag2element<tria_2_tag> : tria_2_elem {};


struct quad_1_tag {};
template <> struct tag2element<quad_1_tag> : quad_1_elem {};
struct quad_2_tag {};
template <> struct tag2element<quad_2_tag> : quad_2_elem {};

struct quad_28_tag {};
template <> struct tag2element<quad_28_tag> : quad_28_elem {};

#endif // LIB_ELEMENT_HPP_INCLUDED
