#include "shapeset.hpp"

template <>
typename constant_shape_set<quad_domain>::shape_t
	const constant_shape_set<quad_domain>::m_shape
	= constant_shape_set<quad_domain>::shape_t::Ones();

template <>
typename constant_shape_set<quad_domain>::dshape_t
	const constant_shape_set<quad_domain>::m_dshape
	= constant_shape_set<quad_domain>::dshape_t::Zero();

template <>
typename constant_shape_set<tria_domain>::shape_t
	const constant_shape_set<tria_domain>::m_shape
	= constant_shape_set<tria_domain>::shape_t::Ones();

template <>
typename constant_shape_set<tria_domain>::dshape_t
	const constant_shape_set<tria_domain>::m_dshape
	= constant_shape_set<tria_domain>::dshape_t::Zero();



quad_2_shape_set::xi_t
	const quad_2_shape_set::m_corners[quad_2_shape_set::num_nodes] = {
		quad_2_shape_set::xi_t(-1.0,-1.0),
		quad_2_shape_set::xi_t( 0.0,-1.0),
		quad_2_shape_set::xi_t(+1.0,-1.0),
		quad_2_shape_set::xi_t(-1.0, 0.0),
		quad_2_shape_set::xi_t( 0.0, 0.0),
		quad_2_shape_set::xi_t(+1.0, 0.0),
		quad_2_shape_set::xi_t(-1.0,+1.0),
		quad_2_shape_set::xi_t( 0.0,+1.0),
		quad_2_shape_set::xi_t(+1.0,+1.0)
};


tria_2_shape_set::xi_t
	const tria_2_shape_set::m_corners[tria_2_shape_set::num_nodes] = {
		tria_2_shape_set::xi_t(0.0, 0.0),
		tria_2_shape_set::xi_t(0.5, 0.0),
		tria_2_shape_set::xi_t(1.0, 0.0),
		tria_2_shape_set::xi_t(0.0, 0.5),
		tria_2_shape_set::xi_t(0.5, 0.5),
		tria_2_shape_set::xi_t(0.0, 1.0)
};

