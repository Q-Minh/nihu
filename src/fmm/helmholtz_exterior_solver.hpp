/**
 * @file helmholtz_exterior_solver.hpp
 * @brief A generic Burton-Miller FMBEM solver for the Helmholtz equation
 * @ingroup fmm_helmholtz
 */

#ifndef HELMHOLTZ_EXTERIOR_SOLVER_HPP_INCLUDED
#define HELMHOLTZ_EXTERIOR_SOLVER_HPP_INCLUDED

#include "core/field.hpp"
#include "core/function_space.hpp"
#include "util/timer.h"
#include "util/type2tag.hpp"

#include "cluster_tree.hpp"
#include "elem_center_iterator.hpp"
#include "fmm_indexed.hpp"
#include "fmm_integrated.hpp"
#include "fmm_matrix.hpp"
#include "fmm_operator_collection.hpp"
#include "fmm_precompute.hpp"
#include "matrix_free.hpp"
#include "preconditioner.hpp"

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
	template <class DivideDerived>
	response_t const &solve(divide_base<DivideDerived> const &divide, size_t far_field_quadrature_order)
	{
		// build the cluster tree
		std::cout << "Building cluster tree ..." << std::endl;
		cluster_tree_t tree(
			create_field_center_iterator(m_trial_space.template field_begin<trial_field_t>()),
			create_field_center_iterator(m_trial_space.template field_end<trial_field_t>()),
			divide.derived());

		// instantiate the fmm object
		std::cout << "Instantiating fmm object ..." << std::endl;
		fmm_t fmm(m_wave_number);
		/// \todo this should be independent of the solver
		fmm.set_accuracy(3.0);

		// initialize tree data
		std::cout << "Initializing tree data ..." << std::endl;
		fmm.init_level_data(tree);
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));
#if PARALLEL
		auto max_num_threads = omp_get_max_threads();
		std::cout << "Expanding to " << max_num_threads << " threads" << std::endl;
		for (size_t i = 0; i < tree.get_n_levels(); ++i)
			fmm.get_level_data(i).set_num_threads(max_num_threads);
#endif

		// create interaction lists
		std::cout << "Computing interaction lists..." << std::endl;
		interaction_lists lists(tree);

		// create functors
		auto int_fctr = create_integrated_functor(test_field_tag_t(), trial_field_tag_t(),
			far_field_quadrature_order, true);

		auto idx_fctr = create_indexed_functor(
			m_test_space.template field_begin<test_field_t>(),
			m_test_space.template field_end<test_field_t>(),
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>(),
			tree);

		auto pre_fctr = create_precompute_functor(tree, lists);

		// Burton-Miller coupling constant
		std::complex<double> alpha(0.0, -1.0 / m_wave_number);

		// integration
		auto I = create_identity_p2p_integral(type2tag<test_field_t>(), type2tag<trial_field_t>());
		auto lhs_collection = create_fmm_operator_collection(
			int_fctr(fmm.template create_p2p<0, 1>())
				+ alpha * int_fctr(fmm.template create_p2p<1, 1>())
				- 0.5 * I,
			int_fctr(fmm.template create_p2m<1>()),
			int_fctr(fmm.template create_p2l<1>()),
			int_fctr(fmm.template create_m2p<0>())
				+ alpha * int_fctr(fmm.template create_m2p<1>()),
			int_fctr(fmm.template create_l2p<0>())
				+ alpha * int_fctr(fmm.template create_l2p<1>()),
			fmm.create_m2m(),
			fmm.create_m2l(),
			fmm.create_l2l()
		);

		auto rhs_collection = create_fmm_operator_collection(
			int_fctr(fmm.template create_p2p<0, 0>())
				+ alpha * int_fctr(fmm.template create_p2p<1, 0>())
				+ (alpha / 2.0) * I,
			int_fctr(fmm.template create_p2m<0>()),
			int_fctr(fmm.template create_p2l<0>())
		);

		// indexing
		auto lhs_cix_collection = lhs_collection.transform(idx_fctr);
		auto rhs_cix_collection = rhs_collection.transform(idx_fctr);

		// precomputation
		std::cout << "Precomputing fmm operators ..." << std::endl;
		auto lhs_pre_collection = lhs_cix_collection.transform(pre_fctr);

		// create rhs matrix object
		std::cout << "Assembling rhs matrix ..." << std::endl;
		auto slp_matrix = create_fmm_matrix(
			pre_fctr(rhs_cix_collection.get(p2p_tag())),
			rhs_cix_collection.get(p2m_tag()),
			rhs_cix_collection.get(p2l_tag()),
			lhs_pre_collection.get(m2p_tag()),
			lhs_pre_collection.get(l2p_tag()),
			lhs_pre_collection.get(m2m_tag()),
			lhs_pre_collection.get(l2l_tag()),
			lhs_pre_collection.get(m2l_tag()),
			tree, lists);

		// compute rhs with fmbem
		std::cout << "Computing rhs ..." << std::endl;
		response_t rhs = slp_matrix * m_excitation;

		// create matrix object
		std::cout << "Assembling lhs matrix ..." << std::endl;
		auto dlp_matrix = create_fmm_matrix(
			lhs_pre_collection,
			tree, lists);

		// compute solution
		std::cout << "Starting iterative solution ..." << std::endl;
		auto M = create_matrix_free(dlp_matrix);

		Eigen::GMRES<decltype(M), Eigen::IdentityPreconditioner > solver(M);
		solver.setTolerance(1e-8);
		solver.set_restart(3000);
		m_response = solver.solve(rhs);
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
} // end of namespace NiHu

#endif /* HELMHOLTZ_EXTERIOR_SOLVER_HPP_INCLUDED */
