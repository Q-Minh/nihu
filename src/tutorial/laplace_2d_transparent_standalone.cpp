#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "../library/lib_element.hpp"

#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

static NiHu::mesh<tmp::vector<NiHu::line_1_elem> > create_mesh(double r, int N)
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

template <class TestSpace, class TrialSpace>
static void tester(TestSpace const &test_space, TrialSpace const &trial_space)
{
	// instantiate integral operators
	auto I_op = NiHu::identity_integral_operator();
	auto G_op = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto H_op = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());
	auto Ht_op = NiHu::create_integral_operator(NiHu::laplace_2d_DLPt_kernel());
	auto D_op = NiHu::create_integral_operator(NiHu::laplace_2d_HSP_kernel());
	
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
			p0(k) += A(s) * -std::log(r)/(2.0 * M_PI);
			q0(k) += A(s) * -(1.0/r) * rdn /(2.0 * M_PI);
		}
	}
	
	std::cout << "Matrix generation" << std::endl;
	
	dMatrix G(N, N), H(N, N), Ht(N, N), D(N, N), I(N,N);
	G.setZero();
	H.setZero();
	Ht.setZero();
	D.setZero();
	I.setZero();

	I << test_space * I_op[trial_space];
	G << test_space * G_op[trial_space];
	H << test_space * H_op[trial_space];
	Ht << test_space * Ht_op[trial_space];
	D << test_space * D_op[trial_space];
	
	// solve with direct Neumann
	std::cout << " Solution with direct Neumann:" << std::endl;
	dMatrix p_direct_neumann = (H - .5 * I).colPivHouseholderQr().solve(G * q0);
	double error_direct_neumann = (p_direct_neumann-p0).norm() / p0.norm();
	std::cout << "Direct Neumann Error: " << error_direct_neumann << std::endl;
	
	// solve with direct Dirichlet
	std::cout << " Solution with direct Dirichlet:" << std::endl;
	dMatrix q_direct_dirichlet = G.colPivHouseholderQr().solve((H - .5 * I) * p0);
	double error_direct_dirichlet = (q_direct_dirichlet-q0).norm() / q0.norm();
	std::cout << "Direct Dirichlet Error: " << error_direct_dirichlet << std::endl;
}

int main(void)
{
	int N = 40;
	double r = 1.2;
	auto mesh = create_mesh(r, N);
	auto const &trial_space = NiHu::constant_view(mesh);
	
	std::cout << "Collocation" << std::endl;
	tester(NiHu::dirac(trial_space), trial_space);
	
	std::cout << "Galerkin" << std::endl;
	tester(trial_space, trial_space);
	
	return 0;
}

