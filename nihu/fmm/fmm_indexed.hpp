/**
 * @file fmm_indexed.hpp 
 * @brief Interface for indexing FMM operators 
 * @ingroup fmm_ops
 */

#ifndef FMM_INDEXED_HPP_INCLUDED
#define FMM_INDEXED_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "p2x_cluster_indexed.hpp"
#include "p2x_indexed.hpp"
#include "x2p_cluster_indexed.hpp"
#include "x2p_indexed.hpp"
#include "p2p_indexed.hpp"
#include "x2x_cluster_indexed.hpp"


#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Op, class FmmTag = typename std::decay<Op>::type::fmm_tag>
class indexed;

template <class Op>
class indexed<Op, m2l_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
	{
		return create_x2x_cluster_indexed(std::forward<Op>(op), tree);
	}
};

template <class Op>
class indexed<Op, l2l_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
	{
		return create_x2x_cluster_indexed(std::forward<Op>(op), tree);
	}
};

template <class Op>
class indexed<Op, m2m_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
	{
		return create_x2x_cluster_indexed(std::forward<Op>(op), tree);
	}
};


template <class Op>
class indexed<Op, m2p_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
	{
		return create_x2p_cluster_indexed(
			create_x2p_indexed(
				std::forward<Op>(op),
				test_begin, test_end
			), tree);
	}
};


template <class Op>
class indexed<Op, l2p_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
	{
		return create_x2p_cluster_indexed(
			create_x2p_indexed(
				std::forward<Op>(op),
				test_begin, test_end
			), tree);
	}
};


template <class Op>
class indexed<Op, p2m_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
	{
		return create_p2x_cluster_indexed(
			create_p2x_indexed(
				std::forward<Op>(op),
				trial_begin, trial_end
			), tree);
	}
};


template <class Op>
class indexed<Op, p2l_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
	{
		return create_p2x_cluster_indexed(
			create_p2x_indexed(
				std::forward<Op>(op),
				trial_begin, trial_end
			), tree);
	}
};


template <class Op>
class indexed<Op, p2p_tag> {
public:
	template <class TestIt, class TrialIt, class ClusterDerived>
	static auto eval(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &)
	{
		return create_p2p_indexed(
			std::forward<Op>(op),
			test_begin, test_end,
			trial_begin, trial_end);
	}
};




template <class Op, class TestIt, class TrialIt, class ClusterDerived>
auto create_indexed(Op &&op, TestIt test_begin, TestIt test_end, TrialIt trial_begin, TrialIt trial_end, cluster_tree<ClusterDerived> const &tree)
{
	return indexed<Op>::eval(std::forward<Op>(op), test_begin, test_end, trial_begin, trial_end, tree);
}


template <class TestIt, class TrialIt, class ClusterDerived>
struct indexed_functor
{
	indexed_functor(TestIt test_begin, TestIt test_end,
		TrialIt trial_begin, TrialIt trial_end,
		cluster_tree<ClusterDerived> const &tree)
		: m_test_begin(test_begin)
		, m_test_end(test_end)
		, m_trial_begin(trial_begin)
		, m_trial_end(trial_end)
		, m_tree(tree)
	{
	}

	template <class Op>
	auto operator()(Op &&op) const
	{
		return create_indexed(std::forward<Op>(op), m_test_begin, m_test_end, m_trial_begin, m_trial_end, m_tree);
	}

	TestIt m_test_begin, m_test_end;
	TrialIt m_trial_begin, m_trial_end;
	cluster_tree<ClusterDerived> const &m_tree;
};

template <class TestIt, class TrialIt, class ClusterDerived>
indexed_functor<TestIt, TrialIt, ClusterDerived>
create_indexed_functor(TestIt test_begin, TestIt test_end,
	TrialIt trial_begin, TrialIt trial_end,
	cluster_tree<ClusterDerived> const &tree)
{
	return indexed_functor<TestIt, TrialIt, ClusterDerived>(
		test_begin, test_end, trial_begin, trial_end,
		tree);
}

}	// end of namespace fmm
}	// end of namespace NiHu

#endif /* FMM_INDEXED_HPP_INCLUDED */
