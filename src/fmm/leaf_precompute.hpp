#ifndef LEAF_PRECOMPUTE_HPP_INCLUDED
#define LEAF_PRECOMPUTE_HPP_INCLUDED

#include <Eigen/SparseCore>

#include "cluster_tree.hpp"
#include "lists.hpp"

#include <vector>
#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator>
class p2x_precompute
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::tree_t tree_t;
	typedef interaction_lists::list_t list_t;
	typedef typename operator_t::result_t result_t;

	static size_t const num_dof_per_src = operator_t::num_dof_per_src;

	p2x_precompute(operator_t const&op, list_t const &list)
		: m_tree(op.get_tree())
		, m_double_idx(m_tree.get_n_clusters(), m_tree.get_n_clusters())
	{
		typedef Eigen::Triplet<size_t, size_t> triplet_t;
		std::vector<triplet_t> triplets;
		size_t idx = 0;
		for (size_t to = 0; to < list.size(); ++to)
		{
			for (auto from : list[to])
			{
				m_container.push_back(op(to, from));
				triplets.push_back(triplet_t(to, from, idx));
				++idx;
			}
		}
		m_double_idx.setFromTriplets(triplets.begin(), triplets.end());
	}

	p2x_precompute(operator_t const &op, std::vector<size_t> const &list)
		: m_tree(op.get_tree())
		, m_single_idx(m_tree.get_n_clusters())
	{
		for (size_t i = 0; i < list.size(); ++i)
		{
			size_t to = list[i];
			m_container.push_back(op(to));
			m_single_idx[to] = i;
		}
	}

	result_t const &operator()(size_t to) const
	{
		return m_container[m_single_idx[to]];
	}

	result_t const &operator()(size_t to, size_t from) const
	{
		return m_container[m_double_idx.coeff(to, from)];
	}

private:
	tree_t const &m_tree;
	std::vector<result_t> m_container;
	Eigen::SparseMatrix<size_t> m_double_idx;
	std::vector<size_t> m_single_idx;
};

template <class Operator>
p2x_precompute<Operator>
create_p2x_precompute(Operator const &op, interaction_lists::list_t const &list)
{
	return p2x_precompute<Operator>(op, list);
}

template <class Operator>
p2x_precompute<Operator>
create_p2x_precompute(Operator const &op, std::vector<size_t> const &list)
{
	return p2x_precompute<Operator>(op, list);
}



template <class Operator>
class x2p_precompute
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::tree_t tree_t;
	typedef interaction_lists::list_t list_t;
	typedef typename operator_t::result_t result_t;

	static size_t const num_dof_per_rec = operator_t::num_dof_per_rec;

	x2p_precompute(operator_t const&op, list_t const &list)
		: m_tree(op.get_tree())
		, m_double_idx(m_tree.get_n_clusters(), m_tree.get_n_clusters())
	{
		typedef Eigen::Triplet<size_t, size_t> triplet_t;
		std::vector<triplet_t> triplets;
		size_t idx = 0;
		for (size_t to = 0; to < list.size(); ++to)
		{
			for (auto from : list[to])
			{
				m_container.push_back(op(to, from));
				triplets.push_back(triplet_t(to, from, idx));
				++idx;
			}
		}
		m_double_idx.setFromTriplets(triplets.begin(), triplets.end());
	}

	x2p_precompute(operator_t const &op, std::vector<size_t> const &list)
		: m_tree(op.get_tree())
		, m_single_idx(m_tree.get_n_clusters())
	{
		for (size_t i = 0; i < list.size(); ++i)
		{
			size_t to = list[i];
			m_container.push_back(op(to));
			m_single_idx[to] = i;
		}
	}

	result_t const &operator()(size_t to) const
	{
		return m_container[m_single_idx[to]];
	}

	result_t const &operator()(size_t to, size_t from) const
	{
		return m_container[m_double_idx.coeff(to, from)];
	}

private:
	tree_t const &m_tree;
	std::vector<result_t> m_container;
	Eigen::SparseMatrix<size_t> m_double_idx;
	std::vector<size_t> m_single_idx;
};

template <class Operator>
x2p_precompute<Operator>
create_x2p_precompute(Operator const &op, interaction_lists::list_t const &list)
{
	return x2p_precompute<Operator>(op, list);
}

template <class Operator>
x2p_precompute<Operator>
create_x2p_precompute(Operator const &op, std::vector<size_t> const &list)
{
	return x2p_precompute<Operator>(op, list);
}


} // end of namespace fmm
} // namespace NiHu

#endif // LEAF_PRECOMPUTE_HPP_INCLUDED
