#include <boost/math/constants/constants.hpp>

#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/laplace_singular_integrals.hpp"

#include <gtest/gtest.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

/// \todo create_circle again
static NiHu::mesh<tmp::vector<NiHu::line_1_elem> > 
create_mesh(double r, int N)
{
	using namespace boost::math::double_constants;
	
	dMatrix nodes(N,2);
	uMatrix elements(N, 3);
	for (int i = 0; i < N; ++i)
	{
		double phi = i*two_pi/N;
		nodes(i,0) = r * cos(phi);
		nodes(i,1) = r * sin(phi);
		elements(i,0) = NiHu::line_1_elem::id;
		elements(i,1) = i % N;
		elements(i,2) = (i+1) % N;
	}
	return NiHu::create_mesh(nodes, elements, NiHu::line_1_tag());
}

template <class TestSpace, class TrialSpace>
static double tester(TestSpace const &test_space, TrialSpace const &trial_space)
{
	using namespace boost::math::double_constants;
	
	// instantiate integral operators
	auto I_op = NiHu::identity_integral_operator();
	auto L_op = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto M_op = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());
	
	unsigned N = test_space.get_num_dofs();
	
	// excitation and response
	dMatrix q0(N, 1), p0(N,1);
	p0.setZero();
	q0.setZero();
	
	dMatrix x0(2, 2);
	x0 <<
		.2, .3,
		.3, .2;
	dMatrix A(1, 2);
	A << 1.0, -1.0;
	
	auto const &mesh = trial_space.get_mesh();
	for (unsigned k = 0; k < N; ++k)
	{
		auto const &elem = mesh.template get_elem<NiHu::line_1_elem>(k);
		auto y = elem.get_center();
		auto ny = elem.get_normal().normalized();
		for (int s = 0; s < x0.cols(); ++s)
		{
			auto rvec = y - x0.col(s);
			double r = rvec.norm();
			double rdn = rvec.dot(ny) / r;
			p0(k) += A(s) * -std::log(r) / two_pi;
			q0(k) += A(s) * -(1.0/r) * rdn / two_pi;
		}
	}
	
	// system matrices
	dMatrix L(N, N), M(N, N), I(N,N);
	L.setZero(); M.setZero(); I.setZero();

	I << test_space * I_op[trial_space];
	L << test_space * L_op[trial_space];
	M << test_space * M_op[trial_space];

	dMatrix p = (M - .5 * I).colPivHouseholderQr().solve(L * q0);

	return (p-p0).norm() / p0.norm();
}

TEST(Laplace2dTransparent, LinearLine_Constant_Collocation)
{
	int N = 26;
	double r = 1.2;
	auto mesh = create_mesh(r, N);
	auto const &trial_space = NiHu::constant_view(mesh);
	
	EXPECT_LE(tester(NiHu::dirac(trial_space), trial_space), 1e-2);
}

TEST(Laplace2dTransparent, LinearLine_Constant_Galerkin)
{
	int N = 26;
	double r = 1.2;
	auto mesh = create_mesh(r, N);
	auto const &trial_space = NiHu::constant_view(mesh);
	
	EXPECT_LE(tester(trial_space, trial_space), 1e-2);
}
