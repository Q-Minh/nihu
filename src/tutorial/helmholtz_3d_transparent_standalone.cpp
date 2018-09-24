#include "core/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"
#include "../library/lib_element.hpp"
#include "../interface/read_off_mesh.hpp"

#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

template <class TestSpace, class TrialSpace>
static double tester(TestSpace const &test_space, TrialSpace const &trial_space)
{
	double wn = 1.;
		
	// instantiate integral operators
	auto I_op = NiHu::identity_integral_operator();
	auto L_op = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(wn));
	auto M_op = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(wn));
	
	unsigned N = test_space.get_num_dofs();
	
	// excitation and response
	cMatrix q0(N, 1), p0(N,1);
	p0.setZero();
	q0.setZero();
	
	dMatrix x0(3, 2);
	x0 <<
		.2, .3,
		.3, .2,
		.1, .2;
	dMatrix A(1, 3);
	A << 1.0, -1.0, -2.0;
	
	auto const &mesh = trial_space.get_mesh();
	for (unsigned k = 0; k < N; ++k)
	{
		auto const &elem = mesh.template get_elem<NiHu::tria_1_elem>(k);
		auto y = elem.get_center();
		auto ny = elem.get_normal().normalized();
		for (int s = 0; s < x0.cols(); ++s)
		{
			auto rvec = y - x0.col(s);
			double r = rvec.norm();
			std::complex<double> ikr(0., wn*r);
			double rdn = rvec.dot(ny) / r;
			p0(k) += A(s) * std::exp(-ikr)/r/(4. * M_PI);
			q0(k) += A(s) * -(std::exp(-ikr)/r/r) * rdn /(4. * M_PI) * (1. + ikr);
		}
	}
	
	// system matrices
	cMatrix L(N, N), M(N, N), I(N,N);
	L.setZero(); M.setZero(); I.setZero();

	I << test_space * I_op[trial_space];
	L << test_space * L_op[trial_space];
	M << test_space * M_op[trial_space];
	
	std::cout << I.block(0, 0, 5, 5) << std::endl << std::endl;
	std::cout << L.block(0, 0, 5, 5) << std::endl << std::endl;
	std::cout << M.block(0, 0, 5, 5) << std::endl << std::endl;

	cMatrix p = (M - .5 * I).colPivHouseholderQr().solve(L * q0);
	
	return (p-p0).norm() / p0.norm();
}

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cerr << "Usage: " << argv[0] << " OFFname" << std::endl;
		return 1;
	}
	
	// read mesh
	auto mesh = NiHu::read_off_mesh(argv[1], NiHu::tria_1_tag());
	auto const &trial_space = NiHu::constant_view(mesh);
	
	std::cout << "Collocation" << std::endl;
	double coll_error = tester(NiHu::dirac(trial_space), trial_space);
	std::cout << coll_error << std::endl;
	
	std::cout << "Galerkin" << std::endl;
	double gal_error = tester(trial_space, trial_space);
	std::cout << gal_error << std::endl;
	
	return 0;
}

