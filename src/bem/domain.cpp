#include "domain.hpp"

template<>
line_domain::xi_t
	const line_domain::m_center =
	line_domain::xi_t::Zero();


template<>
line_domain::corners_t
	const line_domain::m_corners = {
	-line_domain::xi_t::Ones(),
	line_domain::xi_t::Ones()
	};


template<>
tria_domain::xi_t
	const tria_domain::m_center =
	tria_domain::xi_t::Ones()/3.0;


template<>
tria_domain::corners_t const tria_domain::m_corners = {
	tria_domain::xi_t(0.0,0.0),
	tria_domain::xi_t(1.0,0.0),
	tria_domain::xi_t(0.0,1.0)
	};


template<>
quad_domain::xi_t
	const quad_domain::m_center =
	quad_domain::xi_t::Zero();


template<>
quad_domain::corners_t const quad_domain::m_corners = {
	quad_domain::xi_t(-1.0,-1.0),
	quad_domain::xi_t( 1.0,-1.0),
	quad_domain::xi_t( 1.0, 1.0),
	quad_domain::xi_t(-1.0, 1.0)
	};


template<>
brick_domain::corners_t
	const brick_domain::m_corners = {
	brick_domain::xi_t(-1.0,-1.0,-1.0),
	brick_domain::xi_t( 1.0,-1.0,-1.0),
	brick_domain::xi_t( 1.0, 1.0,-1.0),
	brick_domain::xi_t(-1.0, 1.0,-1.0),
	brick_domain::xi_t(-1.0,-1.0, 1.0),
	brick_domain::xi_t( 1.0,-1.0, 1.0),
	brick_domain::xi_t( 1.0, 1.0, 1.0),
	brick_domain::xi_t(-1.0, 1.0, 1.0)
};


template<>
brick_domain::xi_t
	const brick_domain::m_center =
	brick_domain::xi_t::Zero();

