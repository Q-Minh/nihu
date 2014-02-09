#include "lib_domain.hpp"

std::string const line_domain::m_name = "Line domain";
std::string const tria_domain::m_name = "Tria domain";
std::string const quad_domain::m_name = "Quad domain";
std::string const brick_domain::m_name = "Brick domain";

line_domain::corners_t
	const line_domain::m_corners = {
	line_domain::xi_t::Constant(-1.0),
	line_domain::xi_t::Constant(1.0)
	};

line_domain::xi_t const line_domain::m_center = line_domain::xi_t::Zero();

tria_domain::corners_t
    const tria_domain::m_corners = {
	tria_domain::xi_t(0.0,0.0),
	tria_domain::xi_t(1.0,0.0),
	tria_domain::xi_t(0.0,1.0)
	};

tria_domain::xi_t const tria_domain::m_center = tria_domain::xi_t::Constant(1.0/3.0);

quad_domain::corners_t
    const quad_domain::m_corners = {
	quad_domain::xi_t(-1.0,-1.0),
	quad_domain::xi_t( 1.0,-1.0),
	quad_domain::xi_t( 1.0, 1.0),
	quad_domain::xi_t(-1.0, 1.0)
	};

quad_domain::xi_t const quad_domain::m_center = quad_domain::xi_t::Zero();

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

brick_domain::xi_t const brick_domain::m_center = brick_domain::xi_t::Zero();

