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
#include "fmm_indexed.hpp"
#include "fmm_matrix.hpp"
#include "fmm_operator_collection.hpp"
#include "fmm_precompute.hpp"
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

		/// \todo this should be independent of the solver
		fmm.set_accuracy(3.0);


		// build the cluster tree
		std::cout << "Create cluster tree" << std::endl;
		cluster_tree_t tree(
			create_field_center_iterator(m_trial_space.template field_begin<trial_field_t>()),
			create_field_center_iterator(m_trial_space.template field_end<trial_field_t>()),
			divide);
		std::cout << tree << std::endl;

		// initialize tree data
		std::cout << "Initialize level data" << std::endl;


		fmm.init_level_data(tree);
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));

		// create interaction lists
		std::cout << "Compute interaction lists" << std::endl;
		interaction_lists lists(tree);

		// get operators from fmm
		std::cout << "Generate operators" << std::endl;

		auto x2x_collection = create_fmm_operator_collection(
			fmm.create_m2m(),
			fmm.create_m2l(),
			fmm.create_l2l());

		// Apply cluster indexing
		auto idx_fctr = create_indexed_functor(
			m_test_space.template field_begin<test_field_t>(),
			m_test_space.template field_end<test_field_t>(),
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>(),
			tree);
			auto cix_collection = x2x_collection.transform(idx_fctr);

			// Precalculation
			auto pre_fctr = create_precompute_functor(tree, lists);
			auto pre_collection = cix_collection.transform(pre_fctr);

			// integrate operators over fields
			std::cout << "Operator integration" << std::endl;

			auto fmb = create_fmbem(fmm, test_field_tag_t(), trial_field_tag_t(),
				far_field_quadrature_order, true);

			// perform operator arithmetics

			std::complex<double> alpha(0.0, -1.0 / m_wave_number);

			/// ASSEMBLE M2P
			auto im2p_bm = fmb.template create_m2p<0>() + alpha * fmb.template create_m2p<1>();
			auto il2p_bm = fmb.template create_l2p<0>() + alpha * fmb.template create_l2p<1>();

			auto I = create_identity_p2p_integral(type2tag<test_field_t>(), type2tag<trial_field_t>());

			/// ASSEMBLE P2P_BM_0
			auto ip2p_0_bm = fmb.template create_p2p<0, 0>()
				+ alpha * fmb.template create_p2p<1, 0>()
				+ (alpha / 2.0) * I;
			auto ix_p2p_0 = create_p2x_indexed(
				create_x2p_indexed(ip2p_0_bm,
					m_test_space.template field_begin<test_field_t>(),
					m_test_space.template field_end<test_field_t>()),
				m_trial_space.template field_begin<trial_field_t>(),
				m_trial_space.template field_end<trial_field_t>()
			);
			std::cout << "Precomputing P2P_0..." << std::endl;
			auto start = NiHu::wc_time::tic();
			auto p2p_near_0 = p2p_precompute(ix_p2p_0, tree, lists.get_list(lists.P2P));
			std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;


			/// ASSEMBLE P2P_BM_1
			auto ip2p_1_bm = fmb.template create_p2p<0, 1>()
				+ alpha * fmb.template create_p2p<1, 1>()
				- 0.5 * I;
			auto ix_p2p_1 = create_p2x_indexed(
				create_x2p_indexed(ip2p_1_bm,
					m_test_space.template field_begin<test_field_t>(),
					m_test_space.template field_end<test_field_t>()),
				m_trial_space.template field_begin<trial_field_t>(),
				m_trial_space.template field_end<trial_field_t>()
			);
			std::cout << "Precomputing P2P_1..." << std::endl;
			start = NiHu::wc_time::tic();
			auto p2p_near_1 = p2p_precompute(
				ix_p2p_1,
				tree,
				lists.get_list(lists.P2P));
			std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

			// create indexed fmbem operators
			std::cout << "Operator indexing" << std::endl;

			// assign tree to operators for cluster_indexing
			std::cout << "Cluster indexing" << std::endl;

			auto lhs_collection = create_fmm_operator_collection(
				fmb.template create_p2m<1>(),
				fmb.template create_p2l<1>(),
				im2p_bm,
				il2p_bm
			);

			auto rhs_collection = create_fmm_operator_collection(
				fmb.template create_p2m<0>(),
				fmb.template create_p2l<0>()
			);

			auto lhs_cix_collection = lhs_collection.transform(idx_fctr);
			auto rhs_cix_collection = rhs_collection.transform(idx_fctr);

			auto lhs_pre_collection = lhs_cix_collection.transform(pre_fctr);

			// create precomputed fmbem operators

#if PARALLEL
			auto max_num_threads = omp_get_max_threads();
			std::cout << "Expanding to " << max_num_threads << " threads" << std::endl;
			for (size_t i = 0; i < tree.get_n_levels(); ++i)
				fmm.get_level_data(i).set_num_threads(max_num_threads);
#endif

			// create rhs matrix object
			std::cout << "Starting assembling matrices ..." << std::endl;
			auto slp_matrix = create_fmm_matrix(
				p2p_near_0,
				rhs_cix_collection.get(p2m_tag()),
				rhs_cix_collection.get(p2l_tag()),
				lhs_cix_collection.get(m2p_tag()),
				lhs_cix_collection.get(l2p_tag()),
				pre_collection.get(m2m_tag()),
				pre_collection.get(l2l_tag()),
				pre_collection.get(m2l_tag()),
				tree, lists);

			// compute rhs with fmbem
			std::cout << "Computing rhs ..." << std::endl;
			response_t rhs = slp_matrix * m_excitation;
			std::cout << "rhs ready" << std::endl;

			// create matrix object
			auto dlp_matrix = create_fmm_matrix(
				p2p_near_1,
				lhs_pre_collection.get(p2m_tag()),
				lhs_pre_collection.get(p2l_tag()),
				lhs_pre_collection.get(m2p_tag()),
				lhs_pre_collection.get(l2p_tag()),
				pre_collection.get(m2m_tag()),
				pre_collection.get(l2l_tag()),
				pre_collection.get(m2l_tag()),
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
