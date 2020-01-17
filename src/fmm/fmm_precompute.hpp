/**
 * @file fmm_precompute.hpp
 * @brief Operator pre-computation interface 
 * @ingroup fmm_ops
 */ 

#ifndef FMM_PRECOMPUTE_HPP_INCLUDED
#define FMM_PRECOMPUTE_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "leaf_precompute.hpp"
#include "lists.hpp"
#include "p2p_precompute.hpp"
#include "x2x_precompute.hpp"

#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Op, class FmmTag = typename std::decay<Op>::type::fmm_tag>
class precompute;

template <class Op>
class precompute<Op, m2l_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return create_x2x_precompute(std::forward<Op>(op), lists.get_list(lists.M2L));
	}
};


template <class Op>
class precompute<Op, m2m_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return  create_x2x_precompute(std::forward<Op>(op), lists.get_list(lists.M2M));
	}
};


template <class Op>
class precompute<Op, l2l_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return create_x2x_precompute(std::forward<Op>(op), lists.get_list(lists.L2L));
	}
};


template <class Op>
class precompute<Op, p2m_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return create_p2x_precompute(std::forward<Op>(op), tree);
	}
};


template <class Op>
class precompute<Op, p2l_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return create_p2x_precompute(std::forward<Op>(op), tree, lists.get_list(lists.P2L));
	}
};


template <class Op>
class precompute<Op, m2p_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return create_x2p_precompute(std::forward<Op>(op), tree, lists.get_list(lists.M2P));
	}
};


template <class Op>
class precompute<Op, l2p_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return create_x2p_precompute(std::forward<Op>(op), tree);
	}
};


template <class Op>
class precompute<Op, p2p_tag>
{
public:
	template <class ClusterDerived>
	static auto eval(Op &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
	{
		return create_p2p_precompute(std::forward<Op>(op), tree, lists.get_list(lists.P2P));
	}
};


template <class Op, class Tree>
auto create_precompute(Op &&op, Tree const &tree, interaction_lists const &lists)
{
	return precompute<Op>::eval(std::forward<Op>(op), tree, lists);
}

template <class ClusterDerived>
struct precompute_functor
{
	precompute_functor(cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
		: m_tree(tree)
		, m_lists(lists)
	{
	}

	template <class Op>
	auto operator()(Op &&op) const
	{
		return create_precompute(std::forward<Op>(op), m_tree, m_lists);
	}

	cluster_tree<ClusterDerived> const &m_tree;
	interaction_lists const &m_lists;
};

template <class ClusterDerived>
auto create_precompute_functor(cluster_tree<ClusterDerived> const &tree, interaction_lists const &lists)
{
	return precompute_functor<ClusterDerived>(tree, lists);
}

} // end of namespace fmm
} // end of namespace NiHu

#endif /* FMM_PRECOMPUTE_HPP_INCLUDED */
