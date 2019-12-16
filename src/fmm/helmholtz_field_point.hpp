/**
 * @file helmholtz_field_point.hpp 
 * @brief Generic implementation of Helmholtz field point computations
 * @ingroup fmm_helmholtz 
 */

#ifndef HELMHOLTZ_FIELD_POINT_HPP_INCLUDED
#define HELMHOLTZ_FIELD_POINT_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "divide.hpp"
#include "elem_center_iterator.hpp"
#include "fmm_matrix.hpp"
#include "fmm_precompute.hpp"
#include "fmm_indexed.hpp"
#include "fmm_operator_collection.hpp"
#include "fmm_integrated.hpp"
#include "matrix_free.hpp"

#include "core/field.hpp"
#include "core/function_space.hpp"

#include "util/timer.h"
#include "util/type2tag.hpp"

#ifdef NIHU_FMM_PARALLEL
#include <omp.h>
#endif

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
	typedef cluster_tree<cluster_t> cluster_tree_t;

	typedef typename tmp::deref<
		typename tmp::begin<
		typename trial_space_t::field_type_vector_t
		>::type
	>::type trial_field_t;

	enum { num_trial_dofs = trial_field_t::num_dofs };

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
		, m_accuracy(3.0)
		, m_far_field_order(6)
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

	void set_accuracy(double accuracy)
	{
		m_accuracy = accuracy;
	}

	void set_far_field_order(size_t far_field_order)
	{
		m_far_field_order = far_field_order;
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

	template <class DivideDerived>
	response_t const &eval(divide_base<DivideDerived> const &divide)
	{
		// create cluster tree
		std::cout << "Building cluster tree ..." << std::endl;
		cluster_tree_t tree(
			create_field_center_iterator(m_trial_space.template field_begin<trial_field_t>()),
			create_field_center_iterator(m_trial_space.template field_end<trial_field_t>()),
			create_field_center_iterator(m_test_space.template field_begin<test_field_t>()),
			create_field_center_iterator(m_test_space.template field_end<test_field_t>()),
			divide.derived());
		std::cout << tree << std::endl;

		// initialize tree data
		std::cout << "Initializing fmm object and level data ..." << std::endl;
		fmm_t fmm(m_wave_number);
		fmm.set_accuracy(m_accuracy);
		fmm.init_level_data(tree);
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));

		// create interaction lists
		std::cout << "Computing interaction lists ..." << std::endl;
		interaction_lists lists(tree);

		// create functors
		auto int_fctr = create_integrated_functor(test_field_tag_t(), trial_field_tag_t(),
			m_far_field_order, false);

		auto idx_fctr = create_indexed_functor(
			m_test_space.template field_begin<test_field_t>(),
			m_test_space.template field_end<test_field_t>(),
			m_trial_space.template field_begin<trial_field_t>(),
			m_trial_space.template field_end<trial_field_t>(),
			tree);

		auto pre_fctr = create_precompute_functor(tree, lists);

		// operator manipulations
		std::cout << "Precomputing fmm operators ..." << std::endl;
		auto pre_collection = create_fmm_operator_collection(
			src_concatenate(int_fctr(fmm.template create_p2p<0, 0>()), int_fctr(fmm.template create_p2p<0, 1>())),
			src_concatenate(int_fctr(fmm.template create_p2m<0>()), int_fctr(fmm.template create_p2m<1>())),
			src_concatenate(int_fctr(fmm.template create_p2l<0>()), int_fctr(fmm.template create_p2l<1>())),
			int_fctr(fmm.template create_m2p<0>()),
			int_fctr(fmm.template create_l2p<0>()),
			fmm.create_m2m(),
			fmm.create_l2l(),
			fmm.create_m2l()
			).transform(idx_fctr).transform(pre_fctr);

#if NIHU_FMM_PARALLEL
		auto max_num_threads = omp_get_max_threads();
		std::cout << "Expanding to " << max_num_threads << " threads" << std::endl;
		for (size_t i = 0; i < tree.get_n_levels(); ++i)
			fmm.get_level_data(i).set_num_threads(max_num_threads);
#endif

		// create matrix objects
		std::cout << "Starting assembling | SLP | DLP | " << std::endl;
		auto combined_matrix = create_fmm_matrix(pre_collection, tree, lists);
		std::cout << "Combined matrix assembled" << std::endl;

		std::cout << "Computing MVP " << std::endl;
		cvector_t xct(m_psurf.rows() + m_qsurf.rows(), 1);
        for (int i = 0; i < m_psurf.rows()/num_trial_dofs; ++i)
		{
			xct.segment(2 * i * num_trial_dofs, num_trial_dofs) =
				-m_qsurf.segment(i * num_trial_dofs, num_trial_dofs);
			xct.segment((2 * i + 1) * num_trial_dofs, num_trial_dofs) =
				m_psurf.segment(i * num_trial_dofs, num_trial_dofs);
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
	double m_accuracy;
	size_t m_far_field_order;
};

template <class FmmTag, class TestSpace, class TrialSpace>
auto create_helmholtz_field_point(FmmTag, TestSpace const &test_space, TrialSpace const &trial_space)
{
	return helmholtz_field_point<typename NiHu::tag2type<FmmTag>::type, TestSpace, TrialSpace>(
		test_space, trial_space);
}

} // end of namespace fmm
} // end of namespace NiHu

#endif /* HELMHOLTZ_FIELD_POINT_HPP_INCLUDED */
