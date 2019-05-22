#include <fstream>
#include <iostream>
#include <string>

#include "interface/read_off_mesh.hpp"

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
#include "lists.hpp"
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
	NiHu::elastostatics_3d_U_kernel ukernel(nu, mu);
	NiHu::elastostatics_3d_T_kernel tkernel(nu, mu);
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
	NiHu::elastostatics_3d_U_kernel ukernel(nu, mu);

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
		// create fmm
		typedef fmm::black_box_fmm<NiHu::elastostatics_3d_U_kernel> fmm_u_t;
		fmm_u_t fmm(ukernel);

		// create operators
		auto p2p = fmm.create_p2p();
		auto p2m = fmm.create_p2m<0>();
		auto p2l = fmm.create_p2l();
		auto l2p = fmm.create_l2p<0>();
		auto m2p = fmm.create_m2p();
		auto m2m = fmm.create_m2m();
		auto m2l = fmm.create_m2l();
		auto l2l = fmm.create_l2l();

		// integrate operators over the trial field
		size_t quadrature_order = 6;
		auto ip2p = fmm::create_p2p_integral(p2p, fmm::type2tag<test_field_t>(), fmm::type2tag<trial_field_t>(), true);
		auto ip2m = fmm::create_p2x_integral(p2m, quadrature_order, fmm::type2tag<trial_field_t>());
		auto ip2l = fmm::create_p2x_integral(p2l, quadrature_order, fmm::type2tag<trial_field_t>());
		auto il2p = fmm::create_x2p_integral(l2p, quadrature_order, fmm::type2tag<test_field_t>());
		auto im2p = fmm::create_x2p_integral(m2p, quadrature_order, fmm::type2tag<test_field_t>());

		// index integrated operators
		auto ixp2p = fmm::create_p2x_indexed(
			fmm::create_x2p_indexed(
				ip2p,
				test_space.template field_begin<test_field_t>(), test_space.template field_end<test_field_t>()),
			trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());
		auto ixp2m = fmm::create_p2x_indexed(ip2m, trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());
		auto ixp2l = fmm::create_p2x_indexed(ip2l, trial_space.template field_begin<trial_field_t>(), trial_space.template field_end<trial_field_t>());
		auto ixl2p = fmm::create_x2p_indexed(il2p, test_space.template field_begin<test_field_t>(), test_space.template field_end<test_field_t>());
		auto ixm2p = fmm::create_x2p_indexed(im2p, test_space.template field_begin<test_field_t>(), test_space.template field_end<test_field_t>());

		// build cluster_tree
		typedef fmm::cluster_tree<fmm_u_t::cluster_t> tree_t;
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

		// cluster indexing
		auto cixp2m = fmm::create_p2x_cluster_indexed(ixp2m, tree);
		auto cixp2l = fmm::create_p2x_cluster_indexed(ixp2l, tree);
		auto cixl2p = fmm::create_x2p_cluster_indexed(ixl2p, tree);
		auto cixm2p = fmm::create_x2p_cluster_indexed(ixm2p, tree);
		auto cixm2m = fmm::create_x2x_cluster_indexed(m2m, tree);
		auto cixm2l = fmm::create_x2x_cluster_indexed(m2l, tree);
		auto cixl2l = fmm::create_x2x_cluster_indexed(l2l, tree);

		// build interaction lists
		fmm::interaction_lists lists(tree);

		// acceleration
		auto l2p_pre = fmm::create_x2p_precompute(cixl2p, tree.get_leaf_indices());
		std::cout << "l2p done" << std::endl;
		auto p2m_pre = fmm::create_p2x_precompute(cixp2m, tree.get_leaf_indices());
		std::cout << "p2m done" << std::endl;
		auto m2m_pre = fmm::create_x2x_precompute(cixm2m, lists.get_list(lists.M2M));
		std::cout << "m2m done" << std::endl;
		auto l2l_pre = fmm::create_x2x_precompute(cixl2l, lists.get_list(lists.L2L));
		std::cout << "l2l done" << std::endl;
		auto m2l_pre = fmm::create_x2x_precompute(cixm2l, lists.get_list(lists.M2L));
		std::cout << "m2l done" << std::endl;
		auto p2l_pre = fmm::create_p2x_precompute(cixp2l, lists.get_list(lists.P2L));
		std::cout << "p2l done" << std::endl;
		auto m2p_pre = fmm::create_x2p_precompute(cixm2p, lists.get_list(lists.M2P));
		std::cout << "m2p done" << std::endl;
		auto p2p_pre = fmm::p2p_precompute(ixp2p, tree, lists.get_list(lists.Near));
		std::cout << "p2p done" << std::endl;

		auto Umat = fmm::create_fmm_matrix(p2p_pre, p2m_pre, p2l_pre, m2p_pre, l2p_pre, m2m_pre, l2l_pre, m2l_pre, tree, lists);

		rhs = Umat * t;
	}

	{
		NiHu::elastostatics_3d_T_kernel tkernel(nu, mu);
		auto T_op = NiHu::create_integral_operator(tkernel);
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
