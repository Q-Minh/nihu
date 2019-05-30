#include "core/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"
#include "library/lib_mesh.hpp"

#include <gtest/gtest.h>


template <class TestSpace, class TrialSpace>
static double tester(TestSpace const &test_space, TrialSpace const &trial_space)
{
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

	// select wave number
	double k = 1.;
	double r = 1.;

	// instantiate integral operators
	auto I_op = NiHu::identity_integral_operator();
	auto L_op = NiHu::create_integral_operator(NiHu::helmholtz_2d_SLP_kernel<double>(k));
	auto M_op = NiHu::create_integral_operator(NiHu::helmholtz_2d_DLP_kernel<double>(k));
	
	unsigned N = test_space.get_num_dofs();
	
	// excitation
	std::complex<double> q0 = 1.0;
	
	// analytical solution of Neumann problem
	double kr = k * r;
	std::complex<double> p0 = q0
		* NiHu::bessel::H<0,2>(std::complex<double>(kr))
		/ (-k * NiHu::bessel::H<1,2>(std::complex<double>(kr)));

	cMatrix L(N, N), M(N, N), I(N,N);
	L.setZero(); M.setZero(); I.setZero();

	I << test_space * I_op[trial_space];
	L << test_space * L_op[trial_space];
	M << test_space * M_op[trial_space];

	cMatrix q(N,1);
	q.setConstant(q0);

	cMatrix p = (M - .5 * I).colPivHouseholderQr().solve(L * q);

	return (p.array()-p0).matrix().norm() / std::abs(p0) / std::sqrt(N);
}

TEST(Helmholtz2dPulsating, LinearLine_Constant_Collocation)
{
	int N = 20;
	double r = 1.;
	auto mesh = NiHu::create_line_1_circle_mesh(r, N);
	auto const &trial_space = NiHu::constant_view(mesh);
	
	EXPECT_LE(tester(NiHu::dirac(trial_space), trial_space), 1e-2);
}

TEST(Helmholtz2dPulsating, LinearLine_Constant_Galerkin)
{
	int N = 20;
	double r = 1.;
	auto mesh = NiHu::create_line_1_circle_mesh(r, N);
	auto const &trial_space = NiHu::constant_view(mesh);
	
	EXPECT_LE(tester(trial_space, trial_space), 1e-2);
}
