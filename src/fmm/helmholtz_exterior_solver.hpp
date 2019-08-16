/// \file helmholtz_exterior_solver.hpp
/// \brief A generic Burton-Miller FMM solver for the Helmholtz equation

#ifndef HELMHOLTZ_EXTERIOR_SOLVER_HPP_INCLUDED
#define HELMHOLTZ_EXTERIOR_SOLVER_HPP_INCLUDED

#include "core/field.hpp"
#include "core/function_space.hpp"
#include "util/timer.h"
#include "util/type2tag.hpp"

#include "cluster_tree.hpp"
#include "elem_center_iterator.hpp"
#include "fmbem.hpp"
#include "fmm_matrix.hpp"
#include "leaf_precompute.hpp"
#include "matrix_free.hpp"
#include "p2p_precompute.hpp"
#include "p2x_cluster_indexed.hpp"
#include "p2x_indexed.hpp"
#include "preconditioner.hpp"
#include "x2p_cluster_indexed.hpp"
#include "x2p_indexed.hpp"
#include "x2x_cluster_indexed.hpp"
#include "x2x_precompute.hpp"

#include <Eigen/IterativeLinearSolvers>
#include "GMRES.h"

namespace NiHu
{
namespace fmm
{

/// \brief a generic cololocational Burton-Miller solver 
/// \tparam Fmm the fmm type
/// \tparam TrialSpace the type of the trial function space
template <class Fmm, class TrialSpace>
class helmholtz_exterior_solver
{
public:
	typedef Fmm fmm_t;
	typedef TrialSpace trial_space_t;

	typedef dirac_space<trial_space_t> test_space_t;
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> excitation_t;
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> response_t;

	typedef typename fmm_t::cluster_t cluster_t;
	typedef	cluster_tree<cluster_t> cluster_tree_t;

	typedef typename tmp::deref<
		typename tmp::begin<
		typename trial_space_t::field_type_vector_t
		>::type
	>::type trial_field_t;

	typedef type2tag<trial_field_t> trial_field_tag_t;
	typedef dirac_field<trial_field_t> test_field_t;
	typedef type2tag<test_field_t> test_field_tag_t;

	helmholtz_exterior_solver(trial_space_t const &trial_space)
		: m_trial_space(trial_space)
		, m_test_space(dirac(m_trial_space))
		, m_wave_number(0.0)
	{
	}

	void set_wave_number(double k)
	{
		m_wave_number = k;
	}

	double get_wave_number() const
	{
		return m_wave_number;
	}

	response_t const &get_response() const
	{
		return m_response;
	}

	void set_excitation(excitation_t const &xct)
	{
		m_excitation = xct;
	}

	/// \brief solve the BIE
	/// \param [in] divide the cluster division functor
	/// \param [in] far_field_quadrature_order the far field quadrature order
	template <class Divide>
	response_t const &solve(Divide const &divide, size_t far_field_quadrature_order)
	{
		// instantiate the fmm object
		fmm_t fmm(m_wave_number);

		// build the cluster tree
		std::cout << "Create cluster tree" << std::endl;
		cluster_tree_t tree(
			create_field_center_iterator(m_trial_space.template field_begin<trial_field_t>()),
			create_field_center_iterator(m_trial_space.template field_end<trial_field_t>()),
			divide);
		std::cout << tree << std::endl;

		// initialize tree data
		std::cout << "Initialize level data" << std::endl;
		fmm.set_accuracy(3.0);
		fmm.init_level_data(tree);
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));

		// create interaction lists
		std::cout << "Compute interaction lists" << std::endl;
		interaction_lists lists(tree);

		// get operators from fmm
		std::cout << "Generate operators" << std::endl;
		auto m2m = fmm.create_m2m();
		auto l2l = fmm.create_l2l();
		auto m2l = fmm.create_m2l();

		// integrate operators over fields
		std::cout << "Operator integration" << std::endl;

		auto fmb = create_fmbem(fmm, test_field_tag_t(), trial_field_tag_t(),
			far_field_quadrature_order, true);

		auto ip2p_00 = fmb.template create_p2p<0, 0>();
		auto ip2p_01 = fmb.template create_p2p<0, 1>();
		auto ip2p_10 = fmb.template create_p2p<1, 0>();
		auto ip2p_11 = fmb.template create_p2p<1, 1>();

		auto ip2m_0 = fmb.template create_p2m<0>();
		auto ip2l_0 = fmb.template create_p2l<0>();
		auto ip2m_1 = fmb.template create_p2m<1>();
		auto ip2l_1 = fmb.template create_p2l<1>();

		auto im2p_0 = fmb.template create_m2p<0>();
		auto il2p_0 = fmb.template create_l2p<0>();
		auto im2p_1 = fmb.template create_m2p<1>();
		auto il2p_1 = fmb.template create_l2p<1>();

		// perform operator arithmetics
		std::cout << "Operator arithmetics" << std::endl;

		std::complex<double> alpha(0.0, -1.0 / m_wave_number);
		auto im2p_bm = im2p_0 + alpha * im2p_1;
		auto il2p_bm = il2p_0 + alpha * il2p_1;
		auto I = create_identity_p2p_integral(type2tag<test_field_t>(), type2tag<trial_field_t>());
		auto ip2p_0_bm = ip2p_00 + alpha * ip2p_10 + (alpha / 2.0) * I;
		auto ip2p_1_bm = ip2p_01 + alpha * ip2p_11 - 0.5 * I;

		// create indexed fmbem operators
		std::cout << "Operator indexing" << std::endl;

		auto ix_p2m_0 = create_p2x_indexed(ip2m_0,
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>());
		auto ix_p2l_0 = create_p2x_indexed(ip2l_0,
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>());
		auto ix_p2m_1 = create_p2x_indexed(ip2m_1,
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>());
		auto ix_p2l_1 = create_p2x_indexed(ip2l_1,
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>());

		auto ix_m2p_bm = create_x2p_indexed(im2p_bm,
			m_test_space.template field_begin<test_field_t>(),
			m_test_space.template field_end<test_field_t>());
		auto ix_l2p_bm = create_x2p_indexed(il2p_bm,
			m_test_space.template field_begin<test_field_t>(),
			m_test_space.template field_end<test_field_t>());

		auto ix_p2p_0 = create_p2x_indexed(
			create_x2p_indexed(ip2p_0_bm,
				m_test_space.template field_begin<test_field_t>(),
				m_test_space.template field_end<test_field_t>()),
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>()
		);
		auto ix_p2p_1 = create_p2x_indexed(
			create_x2p_indexed(ip2p_1_bm,
				m_test_space.template field_begin<test_field_t>(),
				m_test_space.template field_end<test_field_t>()),
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>()
		);

		// assign tree to operators for cluster_indexing
		std::cout << "Cluster indexing" << std::endl;

		auto cix_m2m = create_x2x_cluster_indexed(m2m, tree);
		auto cix_m2l = create_x2x_cluster_indexed(m2l, tree);
		auto cix_l2l = create_x2x_cluster_indexed(l2l, tree);

		auto cix_p2m_0 = create_p2x_cluster_indexed(ix_p2m_0, tree);
		auto cix_p2m_1 = create_p2x_cluster_indexed(ix_p2m_1, tree);
		auto cix_p2l_0 = create_p2x_cluster_indexed(ix_p2l_0, tree);
		auto cix_p2l_1 = create_p2x_cluster_indexed(ix_p2l_1, tree);

		auto cix_m2p_bm = create_x2p_cluster_indexed(ix_m2p_bm, tree);
		auto cix_l2p_bm = create_x2p_cluster_indexed(ix_l2p_bm, tree);

		// create precomputed fmbem operators
		std::cout << "Precomputing M2M..." << std::endl;
		auto start = NiHu::wc_time::tic();
		auto m2m_pre = create_x2x_precompute(cix_m2m, lists.get_list(lists.M2M));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing L2L..." << std::endl;
		start = NiHu::wc_time::tic();
		auto l2l_pre = create_x2x_precompute(cix_l2l, lists.get_list(lists.L2L));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing M2L..." << std::endl;
		start = NiHu::wc_time::tic();
		auto m2l_pre = create_x2x_precompute(cix_m2l, lists.get_list(lists.M2L));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s\n" << std::endl;

		std::cout << "Precomputing P2M_1..." << std::endl;
		start = NiHu::wc_time::tic();
		auto p2m_1_pre = create_p2x_precompute(cix_p2m_1, tree.get_leaf_src_indices());
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing P2L_1..." << std::endl;
		start = NiHu::wc_time::tic();
		auto p2l_1_pre = create_p2x_precompute(cix_p2l_1, lists.get_list(lists.P2L));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing L2P_BM..." << std::endl;
		start = NiHu::wc_time::tic();
		auto l2p_pre = create_x2p_precompute(cix_l2p_bm, tree.get_leaf_rec_indices());
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing M2P_BM..." << std::endl;
		start = NiHu::wc_time::tic();
		auto m2p_pre = create_x2p_precompute(cix_m2p_bm, lists.get_list(lists.M2P));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing P2P_0..." << std::endl;
		start = NiHu::wc_time::tic();
		auto p2p_near_0 = p2p_precompute(ix_p2p_0, tree, lists.get_list(lists.P2P));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing P2P_1..." << std::endl;
		start = NiHu::wc_time::tic();
		auto p2p_near_1 = p2p_precompute(ix_p2p_1, tree, lists.get_list(lists.P2P));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

#if PARALLEL
		auto max_num_threads = omp_get_max_threads();
		std::cout << "Expanding to " << max_num_threads << " threads" << std::endl;
		for (size_t i = 0; i < tree.get_n_levels(); ++i)
			fmm.get_level_data(i).set_num_threads(max_num_threads);
#endif

		// create rhs matrix object
		std::cout << "Starting assembling matrices ..." << std::endl;
		auto slp_matrix = create_fmm_matrix(
			p2p_near_0, cix_p2m_0, cix_p2l_0, cix_m2p_bm, cix_l2p_bm,
			m2m_pre, l2l_pre, m2l_pre,
			tree, lists);

		// compute rhs with fmbem
		std::cout << "Computing rhs ..." << std::endl;
		response_t rhs = slp_matrix * m_excitation;
		std::cout << "rhs ready" << std::endl;

		// create matrix object
		auto dlp_matrix = create_fmm_matrix(
			p2p_near_1, p2m_1_pre, p2l_1_pre, m2p_pre, l2p_pre,
			m2m_pre, l2l_pre, m2l_pre,
			tree, lists);
		std::cout << "Matrices ready " << std::endl;

		// compute solution
		std::cout << "Starting iterative solution ..." << std::endl;
		auto M = create_matrix_free(dlp_matrix);
		// Eigen::GMRES<decltype(M), diagonal_preconditioner<std::complex<double> > > solver(M);
		Eigen::GMRES<decltype(M), Eigen::IdentityPreconditioner > solver(M);
		solver.setTolerance(1e-8);
		solver.set_restart(3000);
		start = NiHu::wc_time::tic();
		m_response = solver.solve(rhs);
		std::cout << "Ready: " << NiHu::wc_time::toc(start) << " s" << std::endl;
		m_iters = solver.iterations();

		return m_response;
	}

	size_t get_iterations() const
	{
		return m_iters;
	}

private:
	trial_space_t const &m_trial_space;
	test_space_t const &m_test_space;

	double m_wave_number;
	excitation_t m_excitation;
	response_t m_response;
	size_t m_iters;
};

} // end of namespace fmm
} // namespace NiHu

#endif // HELMHOLTZ_EXTERIOR_SOLVER_HPP_INCLUDED
