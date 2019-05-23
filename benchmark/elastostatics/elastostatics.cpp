#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include "interface/read_off_mesh.hpp"

#include <Eigen/IterativeLinearSolvers>


#include "library/lib_element.hpp"
#include "library/elastostatics"

#include "core/function_space.hpp"
#include "core/integral_operator.hpp"
#include "core/weighted_residual.hpp"

#include "black_box_fmm.hpp"
#include "cluster_tree.hpp"
#include "divide.h"
#include "elem_center_iterator.hpp"
#include "fmm_matrix.hpp"
#include "GMRES.h"
#include "lists.hpp"
#include "matrix_free.hpp"
#include "p2p_integral.hpp"
#include "p2x_integral.hpp"
#include "x2p_integral.hpp"
#include "x2p_indexed.hpp"
#include "p2x_indexed.hpp"
#include "p2x_cluster_indexed.hpp"
#include "x2p_cluster_indexed.hpp"
#include "x2x_cluster_indexed.hpp"
#include "p2p_precompute.hpp"
#include "x2x_precompute.hpp"
#include "leaf_precompute.hpp"

typedef std::chrono::high_resolution_clock my_clock_t;

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

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;

void solve(std::string const& meshname, std::string const& excname, std::string const& respname)
{
	auto mesh = NiHu::read_off_mesh(meshname, NiHu::quad_1_tag());
	std::cout << "Mesh read. Number of elements: " << mesh.get_num_elements() << std::endl;

	// create function spaces
	auto const& trial_space = NiHu::constant_view(mesh, NiHu::_3d());
	auto const& test_space = NiHu::dirac(trial_space);
	size_t n = trial_space.get_num_dofs();
	size_t m = test_space.get_num_dofs();
	std::cout << "Function spaces created. Number of DOFs: " << m << " x " << n << std::endl;

	// create kernel
	double const nu = .33;
	double const mu = 1e9;
	elastostatics prob(nu, mu);
	auto ukernel = prob.get_kernel<0, 0>();
	auto tkernel = prob.get_kernel<0, 1>();
	auto U_op = NiHu::create_integral_operator(ukernel);
	auto T_op = NiHu::create_integral_operator(tkernel);

	// read excitation
	dVector t(n, 1), rhs, u;

	{
		std::ifstream ifs(excname);
		if (!ifs)
			throw std::runtime_error("Could not open file " + excname);
		for (size_t i = 0; i < n; ++i)
			ifs >> t(i, 0);
	}

	{
		dMatrix U(n, n);
		U.setZero();
		U << test_space * U_op[trial_space];
		rhs = U * t;
	}

	{
		dMatrix T = -.5 * dMatrix::Identity(n, n);
		T << test_space * T_op[trial_space];
		u = T.colPivHouseholderQr().solve(rhs);
	}

	{
		std::ofstream ofs(respname);
		if (!ofs)
			throw std::runtime_error("Could not open file " + respname);
		for (size_t i = 0; i < n; ++i)
			ofs << u(i, 0) << std::endl;
	}

}

void solve_fmm(std::string const& meshname, std::string const& excname, std::string const& respname)
{
	typedef NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant, NiHu::_3d> trial_field_t;
	typedef NiHu::dirac_field<trial_field_t> test_field_t;

	auto mesh = NiHu::read_off_mesh(meshname, NiHu::quad_1_tag());
	std::cout << "Mesh read. Number of elements: " << mesh.get_num_elements() << std::endl;

	// create function spaces
	auto const& trial_space = NiHu::constant_view(mesh, NiHu::_3d());
	auto const& test_space = NiHu::dirac(trial_space);
	size_t n = trial_space.get_num_dofs();
	size_t m = test_space.get_num_dofs();
	std::cout << "Function spaces created. Number of DOFs: " << m << " x " << n << std::endl;

	// create kernel
	double const nu = .33;
	double const mu = 1e9;
	elastostatics prob(nu, mu);

	// read excitation
	dVector t(n, 1), rhs, u;

	{
		std::ifstream ifs(excname);
		if (!ifs)
			throw std::runtime_error("Could not open file " + excname);
		for (size_t i = 0; i < n; ++i)
			ifs >> t(i, 0);
	}

	// create fmm
	typedef fmm::black_box_fmm<elastostatics> fmm_t;
	fmm_t fmm(prob);

	auto m2m = fmm.create_m2m();
	auto m2l = fmm.create_m2l();
	auto l2l = fmm.create_l2l();

	auto l2p = fmm.create_l2p<0>();
	auto m2p = fmm.create_m2p<0>();

	size_t quadrature_order = 6;
	auto il2p = fmm::create_x2p_integral(l2p, quadrature_order, fmm::type2tag<test_field_t>());
	auto im2p = fmm::create_x2p_integral(m2p, quadrature_order, fmm::type2tag<test_field_t>());

	auto ixl2p = fmm::create_x2p_indexed(il2p, test_space.template field_begin<test_field_t>(), test_space.template field_end<test_field_t>());
	auto ixm2p = fmm::create_x2p_indexed(im2p, test_space.template field_begin<test_field_t>(), test_space.template field_end<test_field_t>());


	// build cluster_tree
	typedef fmm::cluster_tree<fmm_t::cluster_t> tree_t;
	tree_t tree(
		fmm::create_field_center_iterator(trial_space.field_begin<trial_field_t>()),
		fmm::create_field_center_iterator(trial_space.field_end<trial_field_t>()),
		fmm::create_field_center_iterator(test_space.field_begin<test_field_t>()),
		fmm::create_field_center_iterator(test_space.field_end<test_field_t>()),
		fmm::divide_num_nodes(10)
	);

	size_t cheb_order = 5;
	for (size_t i = 0; i < tree.get_n_clusters(); ++i)
		tree[i].set_chebyshev_order(cheb_order);

	auto cixm2m = fmm::create_x2x_cluster_indexed(m2m, tree);
	auto cixm2l = fmm::create_x2x_cluster_indexed(m2l, tree);
	auto cixl2l = fmm::create_x2x_cluster_indexed(l2l, tree);

	auto cixl2p = fmm::create_x2p_cluster_indexed(ixl2p, tree);
	auto cixm2p = fmm::create_x2p_cluster_indexed(ixm2p, tree);

	// build interaction lists
	fmm::interaction_lists lists(tree);

	auto m2m_pre = fmm::create_x2x_precompute(cixm2m, lists.get_list(lists.M2M));
	std::cout << "m2m done" << std::endl;
	auto l2l_pre = fmm::create_x2x_precompute(cixl2l, lists.get_list(lists.L2L));
	std::cout << "l2l done" << std::endl;
	auto m2l_pre = fmm::create_x2x_precompute(cixm2l, lists.get_list(lists.M2L));
	std::cout << "m2l done" << std::endl;

	auto l2p_pre = fmm::create_x2p_precompute(cixl2p, tree.get_leaf_indices());
	std::cout << "l2p done" << std::endl;
	auto m2p_pre = fmm::create_x2p_precompute(cixm2p, lists.get_list(lists.M2P));
	std::cout << "m2p done" << std::endl;

	{
		// create operators
		auto p2p = fmm.create_p2p<0, 0>();
		auto p2m = fmm.create_p2m<0>();
		auto p2l = fmm.create_p2l<0>();

		// integrate operators over the trial field
		auto ip2p = fmm::create_p2p_integral(p2p, fmm::type2tag<test_field_t>(), fmm::type2tag<trial_field_t>(), true);
		auto ip2m = fmm::create_p2x_integral(p2m, quadrature_order, fmm::type2tag<trial_field_t>());
		auto ip2l = fmm::create_p2x_integral(p2l, quadrature_order, fmm::type2tag<trial_field_t>());

		// index integrated operators
		auto ixp2p = fmm::create_p2x_indexed(
			fmm::create_x2p_indexed(
				ip2p,
				test_space.template field_begin<test_field_t>(), test_space.template field_end<test_field_t>()),
			trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());
		auto ixp2m = fmm::create_p2x_indexed(ip2m, trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());
		auto ixp2l = fmm::create_p2x_indexed(ip2l, trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());

		// cluster indexing
		auto cixp2m = fmm::create_p2x_cluster_indexed(ixp2m, tree);
		auto cixp2l = fmm::create_p2x_cluster_indexed(ixp2l, tree);

		// acceleration
		auto p2m_pre = fmm::create_p2x_precompute(cixp2m, tree.get_leaf_indices());
		std::cout << "p2m done" << std::endl;
		auto p2l_pre = fmm::create_p2x_precompute(cixp2l, lists.get_list(lists.P2L));
		std::cout << "p2l done" << std::endl;
		auto p2p_pre = fmm::p2p_precompute(ixp2p, tree, lists.get_list(lists.Near));
		std::cout << "p2p done" << std::endl;

		auto Umat = fmm::create_fmm_matrix(p2p_pre, p2m_pre, p2l_pre, m2p_pre, l2p_pre, m2m_pre, l2l_pre, m2l_pre, tree, lists);

		std::cout << "rhs MVP" << std::endl;
		rhs = Umat * t;
		std::cout << "done" << std::endl;
	}

	{
		// create operators
		auto p2p = fmm.create_p2p<0, 1>();
		auto p2m = fmm.create_p2m<1>();
		auto p2l = fmm.create_p2l<1>();

		// integrate operators over the trial field
		auto ip2p = fmm::create_p2p_integral(p2p, fmm::type2tag<test_field_t>(), fmm::type2tag<trial_field_t>(), true);
		auto ip2m = fmm::create_p2x_integral(p2m, quadrature_order, fmm::type2tag<trial_field_t>());
		auto ip2l = fmm::create_p2x_integral(p2l, quadrature_order, fmm::type2tag<trial_field_t>());

		// index integrated operators
		auto ixp2p = fmm::create_p2x_indexed(
			fmm::create_x2p_indexed(
				ip2p,
				test_space.template field_begin<test_field_t>(), test_space.template field_end<test_field_t>()),
			trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());
		auto ixp2m = fmm::create_p2x_indexed(ip2m, trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());
		auto ixp2l = fmm::create_p2x_indexed(ip2l, trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());

		// cluster indexing
		auto cixp2m = fmm::create_p2x_cluster_indexed(ixp2m, tree);
		auto cixp2l = fmm::create_p2x_cluster_indexed(ixp2l, tree);

		// acceleration
		auto p2m_pre = fmm::create_p2x_precompute(cixp2m, tree.get_leaf_indices());
		std::cout << "p2m done" << std::endl;
		auto p2l_pre = fmm::create_p2x_precompute(cixp2l, lists.get_list(lists.P2L));
		std::cout << "p2l done" << std::endl;
		auto p2p_pre = fmm::p2p_precompute(ixp2p, tree, lists.get_list(lists.Near));
		std::cout << "p2p done" << std::endl;

		auto Tmat = fmm::create_fmm_matrix(p2p_pre, p2m_pre, p2l_pre, m2p_pre, l2p_pre, m2m_pre, l2l_pre, m2l_pre, tree, lists);

		// compute solution
		std::cout << "Starting iterative solution ..." << std::endl;
		auto M = fmm::create_matrix_free(Tmat);
		// Eigen::GMRES<decltype(M), diagonal_preconditioner<std::complex<double> > > solver(M);
		Eigen::GMRES<decltype(M), Eigen::IdentityPreconditioner > solver(M);
		solver.setTolerance(1e-8);
		solver.set_restart(3000);
		auto start = my_clock_t::now();
		u = solver.solve(rhs);
		auto finish = my_clock_t::now();
		auto elapsed = finish - start;
		std::cout << "Ready: " << elapsed.count() << " s" << std::endl;
		auto iters = solver.iterations();
	}

	{
		std::ofstream ofs(respname);
		if (!ofs)
			throw std::runtime_error("Could not open file " + respname);
		for (size_t i = 0; i < n; ++i)
			ofs << u(i, 0) << std::endl;
	}

}

int main(int argc, char const* argv[])
{
	try
	{
		// call: prog.exe meshname xctname respname
		std::string meshname(argv[1]);
		std::string excname(argv[2]);
		std::string respname(argv[3]);

		solve_fmm(meshname, excname, respname);
	}
	catch (std::exception & exc)
	{
		std::cerr << "Standard exception caught: " << exc.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unhandled exception caught" << std::endl;
	}

	return 0;
}
