#include "core/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

NiHu::mesh<tmp::vector<NiHu::line_1_elem> > create_mesh(double r, int N)
{
	dMatrix nodes(N,2);
	uMatrix elements(N, 3);
	for (int i = 0; i < N; ++i)
	{
		double phi = i*2*M_PI/N;
		nodes(i,0) = r * cos(phi);
		nodes(i,1) = r * sin(phi);
		elements(i,0) = NiHu::line_1_elem::id;
		elements(i,1) = i % N;
		elements(i,2) = (i+1) % N;
	}
	return NiHu::create_mesh(nodes, elements, NiHu::line_1_tag());
}

int main()
{
	int N = 20;
	double r = 1.;
	
	auto mesh = create_mesh(r, N);

	double k = 1.;

	auto I_op = NiHu::identity_integral_operator();
	auto L_op = NiHu::create_integral_operator(NiHu::helmholtz_2d_SLP_kernel<double>(k));
	auto M_op = NiHu::create_integral_operator(NiHu::helmholtz_2d_DLP_kernel<double>(k));
	
	std::complex<double> q0 = 1.0;
	
	double kr = k * r;
	std::complex<double> p0 = q0
		* NiHu::bessel::H<0,2>(std::complex<double>(kr))
		/ (-k * NiHu::bessel::H<1,2>(std::complex<double>(kr)));

	auto const &trial_space = NiHu::constant_view(mesh);
	
	{	// Galerkin
		std::cout << "Galerkin with constant elements" << std::endl;

		auto const &test_space = trial_space;
	
		std::cout << "Number of DOF: " << N << std::endl;

		cMatrix L(N, N), M(N, N), I(N,N);
		L.setZero(); M.setZero(); I.setZero();

		I << test_space * I_op[trial_space];
		L << test_space * L_op[trial_space];
		M << test_space * M_op[trial_space];
	
		cMatrix q(N,1);
		q.setConstant(q0);
	
		cMatrix p = (M - .5 * I).colPivHouseholderQr().solve(L * q);
	
		std::cout << "Analytic: " << p0 << std::endl;
		std::cout << "Mean: " << p.mean() << std::endl;
		std::cout << "log10 Error: " << std::log10((p.array()-p0).matrix().norm() / std::abs(p0) / std::sqrt(N)) << std::endl;
	}

	{	// Collocation
		std::cout << "Collocation with constant elements" << std::endl;

		auto const &test_space = NiHu::dirac(trial_space);
	
		std::cout << "Number of DOF: " << N << std::endl;

		cMatrix L(N, N), M(N, N), I(N,N);
		L.setZero(); M.setZero(); I.setZero();

		I << test_space * I_op[trial_space];
		L << test_space * L_op[trial_space];
		M << test_space * M_op[trial_space];
	
		cMatrix q(N,1);
		q.setConstant(q0);
	
		cMatrix p = (M - .5 * I).colPivHouseholderQr().solve(L * q);
	
		std::cout << "Analytic: " << p0 << std::endl;
		std::cout << "Mean: " << p.mean() << std::endl;
		std::cout << "log10 Error: " << std::log10((p.array()-p0).matrix().norm() / std::abs(p0) / std::sqrt(N)) << std::endl;
	}

	return 0;
}
