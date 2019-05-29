#include <gtest/gtest.h>

#include "fmm/black_box_fmm.hpp"
#include "library/elastostatics"

class elastostatics
{
public:
	template <unsigned int Nx, unsigned int Ny>
	struct kernel : public std::conditional<
		Ny == 0,
		NiHu::elastostatics_3d_U_kernel,
		NiHu::elastostatics_3d_T_kernel
	> {};

	template <unsigned int Nx, unsigned int Ny>
	typename kernel<Nx, Ny>::type
		get_kernel() const
	{
		return typename kernel<Nx, Ny>::type(m_nu, m_mu);
	}

	elastostatics(double nu, double mu)
		: m_nu(nu)
		, m_mu(mu)
	{
	}

private:
	double m_nu;
	double m_mu;
};


TEST(elasto_bb_expansion, single_level)
{
	typedef NiHu::fmm::black_box_fmm<elastostatics> fmm_t;
	typedef fmm_t::cluster_t cluster_t;
	typedef cluster_t::bounding_box_t bounding_box_t;
	typedef bounding_box_t::location_t x_t;

	double nu = 1 / 3.;
	double mu = 1;
	elastostatics prob(nu, mu);

	x_t X, Y, XX, YY;
	X << 0.0, 0.0, 0.0;
	Y << 4.0, 0.0, 0.0;
	XX << .5, .5, 0.0;
	YY << 3.5, .5, 0.0;

	size_t chebysev_order = 5;

	bounding_box_t bbX(X, 2.0), bbY(Y, 2.0);
	bounding_box_t bbXX(XX, 1.0), bbYY(YY, 1.0);

	cluster_t CX, CY, CXX, CYY;
	CX.set_bounding_box(bbX);
	CY.set_bounding_box(bbY);
	CXX.set_bounding_box(bbXX);
	CYY.set_bounding_box(bbYY);

	CX.set_chebyshev_order(chebysev_order);
	CY.set_chebyshev_order(chebysev_order);
	CXX.set_chebyshev_order(chebysev_order);
	CYY.set_chebyshev_order(chebysev_order);

	double d = .1;
	x_t _x, _y, nx, ny;
	_x << XX(0) + d, XX(1) - d, 0.0;
	_y << YY(0) - d, YY(1) + d, 0.0;
	nx << 3, 1, 1;
	ny << 1, 0, 0;
	ny = ny.normalized();

	fmm_t f(prob);
	
#define Ny 1

	size_t const Nx = 0;
	

	auto p2p_op = f.create_p2p<Nx, Ny>();
	auto p2m_op = f.create_p2m<Ny>();

	auto p2p_op_0 = f.create_p2p<Nx, 0>();
	auto p2m_op_0 = f.create_p2m<0>();

	auto p2l_op = f.create_p2l<Ny>();
	auto m2p_op = f.create_m2p<Nx>();
	auto l2p_op = f.create_l2p<Nx>();
	auto m2l_op = f.create_m2l();
	auto m2m_op = f.create_m2m();
	auto l2l_op = f.create_l2l();

	decltype(p2p_op)::test_input_t x(_x);

#if Ny == 0
	decltype(p2p_op)::trial_input_t y(_y);
#elif Ny == 1
	decltype(p2p_op)::trial_input_t y(_y, ny);
#endif

	double const eps_req = 1e-10;

	auto p2p = p2p_op(x, y);
	std::cout << "p2p:\n" << p2p << std::endl;

	double eps = 1e-10;

	auto p2p0 = p2p_op_0(x, decltype(p2p_op_0)::trial_input_t(_y - eps * ny));
	auto p2p0d = p2p_op_0(x, decltype(p2p_op_0)::trial_input_t(_y + eps * ny));
	auto p2pnum = ((p2p0d - p2p0) / (2 * eps)).eval();
	std::cout << "p2pnum:\n" << p2pnum << std::endl;
	EXPECT_LE((p2pnum - p2p).norm() / p2p.norm(), 1e-10);

	// test p2m and m2p (multipole expansion)

	auto p2m = p2m_op(CY, y);

	auto p2m0 = p2m_op_0(CY, decltype(p2m_op_0)::trial_input_t(_y));
	auto p2m0d = p2m_op_0(CY, decltype(p2m_op_0)::trial_input_t(_y + eps * ny));
	auto p2mnum = ((p2m0d - p2m0) / eps).eval();
	EXPECT_LE((p2mnum - p2m).norm() / p2m.norm(), 1e-10);

	auto m2p = m2p_op(x, CY);

	std::cout << "p2m2p:\n" << m2p * p2m << std::endl;

	double eps_p2m2p = (p2p - m2p * p2m).norm() / p2p.norm();
	EXPECT_LE(eps_p2m2p, eps_req);

	// test p2l and l2p (local expansion)
	auto p2l = p2l_op(CX, y);
	auto l2p = l2p_op(x, CX);

	double eps_p2l2p = (p2p - l2p * p2l).norm() / p2p.norm();
	EXPECT_LE(eps_p2l2p, eps_req);

	// test translation
	auto m2l = m2l_op(CX, CY);
	double eps_p2m2l2p = (p2p - l2p * m2l * p2m).norm() / p2p.norm();
	EXPECT_LE(eps_p2m2l2p, eps_req);
	
#if 0
	// test shiftup
	p2m = p2m_op(CYY, y);
	auto m2m = m2m_op(CY, CYY);

	double eps_p2m2m2l2p = std::abs(p2p - l2p * m2l * m2m * p2m).norm() / p2p.norm();
	EXPECT_LE(eps_p2m2m2l2p, eps_req);

	// test shiftdown
	l2p = l2p_op(x, CXX);
	auto l2l = l2l_op(CXX, CX);
	m2l = m2l_op(CX, CY);

	double eps_p2m2m2l2l2p = (p2p - l2p * l2l * m2l * m2m * p2m).norm() / p2p.norm();
	EXPECT_LE(eps_p2m2m2l2l2p, eps_req);
#endif
}
