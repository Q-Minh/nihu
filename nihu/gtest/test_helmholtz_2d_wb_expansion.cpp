#include <gtest/gtest.h>

#include "fmm/helmholtz_2d_wb_fmm.hpp"

TEST(helmholtz_2d_wb_expansion, single_level)
{
	typedef double wave_number_t;
	typedef NiHu::fmm::helmholtz_2d_wb_fmm<wave_number_t> fmm_t;
	typedef fmm_t::cluster_t cluster_t;

	typedef NiHu::fmm::bounding_box<2, double> bounding_box_t;
	typedef bounding_box_t::location_t x_t;

	x_t X, Y, XX, YY;
	X << 0.0, 0.0;
	Y << 4.0, 0.0;
	XX << .5, .5;
	YY << 3.5, .5;

	double d = .1;
	x_t _x, _y, nx, ny;
	_x << XX(0) + d, XX(1) - d;
	_y << YY(0) - d, YY(1) + d;
	nx << 3, 1;
	ny << -1, 2;

	bounding_box_t bbX(X, 2.0), bbY(Y, 2.0);
	bounding_box_t bbXX(XX, 1.0), bbYY(YY, 1.0);


	NiHu::fmm::helmholtz_2d_wb_level_data d0, d1;
	wave_number_t k = 15;
	using namespace boost::math::double_constants;
	double lambda = two_pi / std::real(k);
	d0.init(2.0 / lambda);
	d1.init(1.0 / lambda);

	cluster_t CX, CY, CXX, CYY;
	CX.set_bounding_box(bbX);
	CY.set_bounding_box(bbY);
	CX.set_p_level_data(&d0);
	CY.set_p_level_data(&d0);

	CXX.set_bounding_box(bbXX);
	CYY.set_bounding_box(bbYY);
	CXX.set_p_level_data(&d1);
	CYY.set_p_level_data(&d1);

	fmm_t f(k);
	
	size_t const Ny = 1;
	size_t const Nx = 1;
	
	auto p2p_op = f.create_p2p<Nx, Ny>();
	auto p2m_op = f.create_p2m<Ny>();
	auto p2l_op = f.create_p2l<Ny>();
	auto m2p_op = f.create_m2p<Nx>();
	auto l2p_op = f.create_l2p<Nx>();
	auto m2l_op = f.create_m2l();
	auto m2m_op = f.create_m2m();
	auto l2l_op = f.create_l2l();

	decltype(p2p_op)::test_input_t x(_x, nx);
	decltype(p2p_op)::trial_input_t y(_y, ny);

	auto p2p = p2p_op(x, y);

	auto p2m = p2m_op(CY, y);
	auto p2l = p2l_op(CX, y);
	auto m2p = m2p_op(x, CY);
	auto l2p = l2p_op(x, CX);
	auto m2l = m2l_op(CX, CY);

	double const eps_req = 1e-6;

	// test p2m and m2p (multipole expansion)
	double eps_p2m2p = std::abs(p2p - (m2p * p2m)(0)) / std::abs(p2p);
	EXPECT_LE(eps_p2m2p, eps_req);

	// test p2l and l2p (local expansion)
	double eps_p2l2p = std::abs(p2p - (l2p * p2l)(0)) / std::abs(p2p);
	EXPECT_LE(eps_p2l2p, eps_req);

	// test translation
	double eps_p2m2l2p = std::abs(p2p - (l2p * (m2l * p2m))(0)) / std::abs(p2p);
	EXPECT_LE(eps_p2m2l2p, eps_req);

	// test shiftup and shiftdown
	p2m = p2m_op(CYY, y);
	l2p = l2p_op(x, CXX);
	auto m2m = m2m_op(CY, CYY);
	auto l2l = l2l_op(CXX, CX);
	m2l = m2l_op(CX, CY);

	double eps_p2m2m2l2l2p = std::abs(p2p - (l2p * (l2l * (m2l * (m2m * p2m))))(0)) / std::abs(p2p);
	EXPECT_LE(eps_p2m2m2l2l2p, eps_req);
}
