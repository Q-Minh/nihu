#ifndef HELMHOLTZ_FIELD_POINT_HPP_INCLUDED
#define HELMHOLTZ_FIELD_POINT_HPP_INCLUDED

#include "core/field.hpp"
#include "core/function_space.hpp"

#include "cluster_tree.hpp"
#include "divide.h"
#include "elem_center_iterator.hpp"
#include "fmm_matrix.hpp"
#include "fmm_precompute.hpp"
#include "fmm_indexed.hpp"
#include "fmm_operator_collection.hpp"
#include "fmm_integrated.hpp"
#include "matrix_free.hpp"

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
	typedef cluster_tree<cluster_t> cluster_tree_t;

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
		// create cluster tree
		std::cout << "Building cluster tree ..." << std::endl;
		cluster_tree_t tree(
			create_field_center_iterator(m_trial_space.template field_begin<trial_field_t>()),
			create_field_center_iterator(m_trial_space.template field_end<trial_field_t>()),
			create_field_center_iterator(m_test_space.template field_begin<test_field_t>()),
			create_field_center_iterator(m_test_space.template field_end<test_field_t>()),
			divide);
		std::cout << tree << std::endl;

		// initialize tree data
		std::cout << "Initializing fmm object and level data ..." << std::endl;
		fmm_t fmm(m_wave_number);
		fmm.set_accuracy(3.0);
		fmm.init_level_data(tree);
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));

		// create interaction lists
		std::cout << "Computing interaction lists ..." << std::endl;
		interaction_lists lists(tree);

		// create functors
		auto int_fctr = create_integrated_functor(test_field_tag_t(), trial_field_tag_t(),
			far_field_quadrature_order, false);

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

#if PARALLEL
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
