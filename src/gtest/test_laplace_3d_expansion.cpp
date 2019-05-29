#include <gtest/gtest.h>

#include "fmm/laplace_3d_fmm.hpp"

TEST(laplace_3d_expansion, single_level)
{
	typedef NiHu::fmm::laplace_3d_fmm fmm_t;
	typedef fmm_t::cluster_t cluster_t;

	typedef NiHu::fmm::bounding_box<3, double> bounding_box_t;
	typedef bounding_box_t::location_t x_t;

	x_t X, Y, XX, YY;
	X << 0.0, 0.0, 0.0;
	Y << 4.0, 0.0, 0.0;
	XX << .5, .5, 0.0;
	YY << 3.5, .5, 0.0;

	double d = .1;
	x_t _x, _y, nx, ny;
	_x << XX(0) + d, XX(1) - d, 0.0;
	_y << YY(0) - d, YY(1) + d, 0.0;
	nx << 3, 1, 1;
	ny << -1, 2, 1;

	bounding_box_t bbX(X, 2.0), bbY(Y, 2.0);
	bounding_box_t bbXX(XX, 1.0), bbYY(YY, 1.0);

	size_t expansion_length = 3;

	cluster_t CX, CY, CXX, CYY;
	CX.set_bounding_box(bbX);
	CXX.set_bounding_box(bbXX);
	CX.set_expansion_length(expansion_length);
	CXX.set_expansion_length(expansion_length);
	CY.set_bounding_box(bbY);
	CYY.set_bounding_box(bbYY);
	CY.set_expansion_length(expansion_length);
	CYY.set_expansion_length(expansion_length);

	fmm_t f;
	
	size_t const Ny = 0;
	size_t const Nx = 0;
	
	auto p2p_op = f.create_p2p<Nx, Ny>();
	auto p2m_op = f.create_p2m<Ny>();
	auto p2l_op = f.create_p2l<Ny>();
	auto m2p_op = f.create_m2p<Nx>();
	auto l2p_op = f.create_l2p<Nx>();
	auto m2l_op = f.create_m2l();
	auto m2m_op = f.create_m2m();
	auto l2l_op = f.create_l2l();

	decltype(p2p_op)::test_input_t x(_x);
	decltype(p2p_op)::trial_input_t y(_y);

	auto p2p = p2p_op(x, y);

	auto p2m = p2m_op(CY, y);
	auto p2l = p2l_op(CX, y);
	auto m2p = m2p_op(x, CY);
	auto l2p = l2p_op(x, CX);
	auto m2l = m2l_op(CX, CY);

	double const eps_req = 1e-10;

	// test p2m and m2p (multipole expansion)
	double eps_p2m2p = std::abs(p2p - (m2p * p2m)(0)) / std::abs(p2p);
	EXPECT_LE(eps_p2m2p, eps_req);

	// test p2l and l2p (local expansion)
	double eps_p2l2p = std::abs(p2p - (l2p * p2l)(0)) / std::abs(p2p);
	EXPECT_LE(eps_p2l2p, eps_req);

	// test translation
	double eps_p2m2l2p = std::abs(p2p - (l2p * m2l * p2m)(0).real()) / std::abs(p2p);
	EXPECT_LE(eps_p2m2l2p, eps_req);
	
	// test shiftup
	p2m = p2m_op(CYY, y);
	auto m2m = m2m_op(CY, CYY);

	std::cout << "m2m: " << std::endl;
	std::cout << m2m.cwiseAbs() << std::endl;

	double eps_p2m2m2l2p = std::abs(p2p - (l2p * m2l * m2m * p2m)(0).real()) / std::abs(p2p);
	EXPECT_LE(eps_p2m2m2l2p, eps_req);

	// test shiftdown
	l2p = l2p_op(x, CXX);
	auto l2l = l2l_op(CXX, CX);
	m2l = m2l_op(CX, CY);

	double eps_p2m2m2l2l2p = std::abs(p2p - (l2p * l2l * m2l * m2m * p2m)(0).real()) / std::abs(p2p);
	EXPECT_LE(eps_p2m2m2l2l2p, eps_req);
}
