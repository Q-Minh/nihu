#ifndef HELMHOLTZ_FIELD_POINT_HPP_INCLUDED
#define HELMHOLTZ_FIELD_POINT_HPP_INCLUDED

#include "core/field.hpp"
#include "core/function_space.hpp"

#include "cluster_tree.hpp"
#include "divide.h"
#include "elem_center_iterator.hpp"
#include "fmm_matrix.hpp"
#include "matrix_free.hpp"
#include "p2p_integral.hpp"
#include "p2p_precompute.hpp"
#include "p2x_cluster_indexed.hpp"
#include "p2x_indexed.hpp"
#include "p2x_integral.hpp"
#include "x2p_cluster_indexed.hpp"
#include "x2p_indexed.hpp"
#include "x2p_integral.hpp"
#include "x2x_cluster_indexed.hpp"
#include "x2x_precompute.hpp"

#include "util/timer.h"


namespace NiHu
{
namespace fmm
{

template <class Fmm, class TestSpace, class TrialSpace>
class helmholtz_field_point
{
public:
	typedef Fmm fmm_t;
	typedef TestSpace test_space_t;
	typedef TrialSpace trial_space_t;

	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
	typedef cvector_t excitation_t;
	typedef cvector_t response_t;

	typedef typename fmm_t::cluster_t cluster_t;
	typedef fmm::cluster_tree<cluster_t> cluster_tree_t;

	typedef typename tmp::deref<
		typename tmp::begin<
		typename trial_space_t::field_type_vector_t
		>::type
	>::type trial_field_t;

	typedef typename tmp::deref<
		typename tmp::begin<
		typename test_space_t::field_type_vector_t
		>::type
	>::type test_field_t;

	typedef type2tag<trial_field_t> trial_field_tag_t;
	typedef type2tag<test_field_t> test_field_tag_t;

	helmholtz_field_point(test_space_t const &test_space,
		trial_space_t const &trial_space)
		: m_test_space(test_space)
		, m_trial_space(trial_space)
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

	void set_psurf(excitation_t const &psurf)
	{
		m_psurf = psurf;
	}

	void set_qsurf(excitation_t const &qsurf)
	{
		m_qsurf = qsurf;
	}

	template <class Divide>
	response_t const &eval(Divide const &divide, size_t far_field_quadrature_order)
	{
		fmm_t fmm(m_wave_number);

		// create cluster tree
		std::cout << "Create cluster tree" << std::endl;
		cluster_tree_t tree(
			fmm::create_field_center_iterator(m_trial_space.template field_begin<trial_field_t>()),
			fmm::create_field_center_iterator(m_trial_space.template field_end<trial_field_t>()),
			fmm::create_field_center_iterator(m_test_space.template field_begin<test_field_t>()),
			fmm::create_field_center_iterator(m_test_space.template field_end<test_field_t>()),
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
		fmm::interaction_lists lists(tree);

		// get operators from fmm
		std::cout << "Generate operators" << std::endl;
		auto m2m = fmm.create_m2m();
		auto l2l = fmm.create_l2l();
		auto m2l = fmm.create_m2l();

		// integrate operators over fields
		std::cout << "Operator integration" << std::endl;

		auto fmb = create_fmbem(fmm, test_field_tag_t(), trial_field_tag_t(),
			far_field_quadrature_order, false);

		auto ip2p_00 = fmb.template create_p2p<0, 0>();
		auto ip2p_01 = fmb.template create_p2p<0, 1>();

		auto ip2m_0 = fmb.template create_p2m<0>();
		auto ip2l_0 = fmb.template create_p2l<0>();
		auto ip2m_1 = fmb.template create_p2m<1>();
		auto ip2l_1 = fmb.template create_p2l<1>();

		auto im2p = fmb.template create_m2p<0>();
		auto il2p = fmb.template create_l2p<0>();

		// concatenate operators
		std::cout << "Operator concatenation" << std::endl;

		auto ip2p = src_concatenate(ip2p_00, ip2p_01);
		auto ip2m = src_concatenate(ip2m_0, ip2m_1);
		auto ip2l = src_concatenate(ip2l_0, ip2l_1);

		// create indexed fmbem operators
		std::cout << "Operator indexing" << std::endl;

		auto ix_p2m = fmm::create_p2x_indexed(ip2m,
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>());
		auto ix_p2l = fmm::create_p2x_indexed(ip2l,
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>());

		auto ix_m2p_0 = fmm::create_x2p_indexed(im2p,
			m_test_space.template field_begin<test_field_t>(),
			m_test_space.template field_end<test_field_t>());
		auto ix_l2p_0 = fmm::create_x2p_indexed(il2p,
			m_test_space.template field_begin<test_field_t>(),
			m_test_space.template field_end<test_field_t>());

		auto ix_p2p = fmm::create_p2x_indexed(
			fmm::create_x2p_indexed(ip2p,
				m_test_space.template field_begin<test_field_t>(),
				m_test_space.template field_end<test_field_t>()),
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>()
		);

		// assign tree to operators for cluster_indexing
		std::cout << "Cluster indexing" << std::endl;

		auto cix_m2m = create_x2x_cluster_indexed(fmm.create_m2m(), tree);
		auto cix_m2l = create_x2x_cluster_indexed(fmm.create_m2l(), tree);
		auto cix_l2l = create_x2x_cluster_indexed(fmm.create_l2l(), tree);

		auto cix_p2m = create_p2x_cluster_indexed(ix_p2m, tree);
		auto cix_p2l = create_p2x_cluster_indexed(ix_p2l, tree);

		auto cix_m2p_0 = create_x2p_cluster_indexed(ix_m2p_0, tree);
		auto cix_l2p_0 = create_x2p_cluster_indexed(ix_l2p_0, tree);

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
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;

		std::cout << "Precomputing P2P..." << std::endl;
		start = NiHu::wc_time::tic();
		auto p2p_near = p2p_precompute(ix_p2p, tree, lists.get_list(lists.P2P));
		std::cout << "Ready, Elapsed time: " << NiHu::wc_time::toc(start) << " s" << std::endl;


#if PARALLEL
		auto max_num_threads = omp_get_max_threads();
		std::cout << "Expanding to " << max_num_threads << " threads" << std::endl;
		for (size_t i = 0; i < tree.get_n_levels(); ++i)
			fmm.get_level_data(i).set_num_threads(max_num_threads);
#endif

		// create matrix objects
		std::cout << "Starting assembling | SLP | DLP | " << std::endl;
		auto combined_matrix = fmm::create_fmm_matrix(
			p2p_near, cix_p2m, cix_p2l, cix_m2p_0, cix_l2p_0,
			m2m_pre, l2l_pre, m2l_pre,
			tree, lists);
		std::cout << "Combined matrix assembled" << std::endl;

		std::cout << "Computing MVP " << std::endl;
		cvector_t xct(m_psurf.rows() + m_qsurf.rows(), 1);
		for (int i = 0; i < m_psurf.rows(); ++i)
		{
			xct(2 * i, 0) = -m_qsurf(i, 0);
			xct(2 * i + 1, 0) = m_psurf(i, 0);
		}

		auto wc_t0 = NiHu::wc_time::tic();
		auto cpu_t0 = NiHu::cpu_time::tic();
		m_response = combined_matrix * xct;
		std::cout << "MVP wall clock time: " << NiHu::wc_time::toc(wc_t0) << " s" << std::endl;
		std::cout << "MVP CPU time: " << NiHu::cpu_time::toc(cpu_t0) << " s" << std::endl;

		combined_matrix.get_timer().print(std::cout);

		return m_response;
	}

private:
	test_space_t const &m_test_space;
	trial_space_t const &m_trial_space;

	double m_wave_number;
	excitation_t m_psurf;
	excitation_t m_qsurf;
	response_t m_response;
};

} // end of namespace fmm
} // namespace NiHu

#endif // HELMHOLTZ_2D_FIELD_POINT_HPP_INCLUDED