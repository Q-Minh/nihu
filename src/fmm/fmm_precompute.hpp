#ifndef FMM_PRECOMPUTE_HPP_INCLUDED
#define FMM_PRECOMPUTE_HPP_INCLUDED

#include "leaf_precompute.hpp"
#include "x2x_precompute.hpp"
#include "lists.hpp"

namespace NiHu
{
namespace fmm
{

template <class Op, class FmmTag = typename Op::fmm_tag>
class precompute;

template <class Op>
class precompute<Op, m2l_tag>
{
public:
	typedef x2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op const &op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(op, lists.get_list(lists.M2L));
	}
};


template <class Op>
class precompute<Op, m2m_tag>
{
public:
	typedef x2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op const &op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(op, lists.get_list(lists.M2M));
	}
};


template <class Op>
class precompute<Op, l2l_tag>
{
public:
	typedef x2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op const &op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(op, lists.get_list(lists.L2L));
	}
};


template <class Op>
class precompute<Op, p2m_tag>
{
public:
	typedef p2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op const &op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(op, tree.get_leaf_src_indices());
	}
};



template <class Op>
class precompute<Op, p2l_tag>
{
public:
	typedef p2x_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op const &op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(op, lists.get_list(lists.P2L));
	}
};


template <class Op>
class precompute<Op, m2p_tag>
{
public:
	typedef x2p_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op const &op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(op, lists.get_list(lists.M2P));
	}
};



template <class Op>
class precompute<Op, l2p_tag>
{
public:
	typedef x2p_precompute<Op> return_type;

	template <class Tree>
	static return_type eval(Op const &op, Tree const &tree, interaction_lists const &lists)
	{
		return return_type(op, tree.get_leaf_rec_indices());
	}
};



template <class Op, class Tree>
typename precompute<Op>::return_type
create_precompute(Op const &op, Tree const &tree, interaction_lists const &lists)
{
	return precompute<Op>::eval(op, tree, lists);
}

}
}

#endif
