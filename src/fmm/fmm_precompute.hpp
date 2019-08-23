#ifndef FMM_PRECOMPUTE_HPP_INCLUDED
#define FMM_PRECOMPUTE_HPP_INCLUDED

#include "leaf_precompute.hpp"
#include "x2x_precompute.hpp"
#include "lists.hpp"

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
	typedef x2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op &&op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(std::forward<Op>(op), lists.get_list(lists.M2L));
	}
};


template <class Op>
class precompute<Op, m2m_tag>
{
public:
	typedef x2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op &&op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(std::forward<Op>(op), lists.get_list(lists.M2M));
	}
};


template <class Op>
class precompute<Op, l2l_tag>
{
public:
	typedef x2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op &&op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(std::forward<Op>(op), lists.get_list(lists.L2L));
	}
};


template <class Op>
class precompute<Op, p2m_tag>
{
public:
	typedef p2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op &&op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(std::forward<Op>(op), tree.get_leaf_src_indices());
	}
};


template <class Op>
class precompute<Op, p2l_tag>
{
public:
	typedef p2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op &&op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(std::forward<Op>(op), lists.get_list(lists.P2L));
	}
};


template <class Op>
class precompute<Op, m2p_tag>
{
public:
	typedef x2p_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op &&op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(std::forward<Op>(op), lists.get_list(lists.M2P));
	}
};


template <class Op>
class precompute<Op, l2p_tag>
{
public:
	typedef x2p_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op &&op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(std::forward<Op>(op), tree.get_leaf_rec_indices());
	}
};


template <class Op, class Tree>
typename precompute<Op>::return_type
create_precompute(Op &&op, Tree const &tree, interaction_lists const &lists)
{
	return precompute<Op>::eval(std::forward<Op>(op), tree, lists);
}

template <class Tree, class Lists>
struct precompute_functor
{
	precompute_functor(Tree const &tree, Lists const &lists)
		: m_tree(tree)
		, m_lists(lists)
	{
	}

	template <class Op>
	auto operator()(Op &&op) const
	{
		return create_precompute(std::forward<Op>(op), m_tree, m_lists);
	}

	Tree const &m_tree;
	Lists const &m_lists;
};

template <class Tree, class Lists>
precompute_functor<Tree, Lists>
create_precompute_functor(Tree const &tree, Lists const &lists)
{
	return precompute_functor<Tree, Lists>(tree, lists);
}

}
}

#endif
