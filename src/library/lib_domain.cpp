#include "lib_domain.hpp"

namespace domain_traits
{
	template <>
	const std::string name<line_domain>::value = "Line domain";
	template <>
	std::string const name<tria_domain>::value = "Tria domain";
	template <>
	std::string const name<quad_domain>::value = "Quad domain";
	template <>
	std::string const name<brick_domain>::value = "Brick domain";
}

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

