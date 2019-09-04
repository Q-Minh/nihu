#ifndef LEAF_PRECOMPUTE_HPP_INCLUDED
#define LEAF_PRECOMPUTE_HPP_INCLUDED

#include <Eigen/SparseCore>

#include "cluster_tree.hpp"
#include "fmm_operator.hpp"
#include "lists.hpp"

#include <vector>
#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Result, class FmmTag>
class p2x_precompute;


template <class Result>
class p2x_precompute<Result, p2m_tag>
	: public fmm_operator<p2m_tag>
{
public:
	typedef Result result_t;

	template <class Operator, class ClusterDerived>
	p2x_precompute(Operator const &op,
		cluster_tree<ClusterDerived> const &tree)
		: m_single_idx(tree.get_n_clusters())
	{
		auto const &src_indices = tree.get_leaf_src_indices();
		for (size_t i = 0; i < src_indices.size(); ++i)
		{
			size_t to = src_indices[i];
			m_container.push_back(op(to));
			m_single_idx[to] = i;
		}
	}

	result_t const &operator()(size_t to) const
	{
		return m_container[m_single_idx[to]];
	}

private:
	std::vector<result_t> m_container;
	std::vector<size_t> m_single_idx;
};


template <class Result>
class p2x_precompute<Result, p2l_tag>
	: public fmm_operator<p2l_tag>
{
public:
	typedef Result result_t;
	typedef interaction_lists::list_t list_t;

	template <class Operator, class ClusterDerived>
	p2x_precompute(Operator const &op,
		cluster_tree<ClusterDerived> const &tree,
		list_t const &list)
		: m_double_idx(tree.get_n_clusters(), tree.get_n_clusters())
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

	result_t const &operator()(size_t to, size_t from) const
	{
		return m_container[m_double_idx.coeff(to, from)];
	}

private:
	std::vector<result_t> m_container;
	Eigen::SparseMatrix<size_t> m_double_idx;
};



template <class Operator>
auto create_p2x_precompute(
	Operator const &op,
	cluster_tree<typename std::decay<Operator>::type::cluster_t> const &tree,
	interaction_lists::list_t const &list)
{
	typedef typename std::decay<Operator>::type operator_t;
	return p2x_precompute<
		operator_t::result_t,
		operator_t::fmm_tag
	>(op, tree, list);
}

template <class Operator>
auto create_p2x_precompute(Operator const &op,
	cluster_tree<typename std::decay<Operator>::type::cluster_t> const &tree)
{
	typedef typename std::decay<Operator>::type operator_t;
	return p2x_precompute<
		operator_t::result_t,
		operator_t::fmm_tag
	>(op, tree);
}


template <class Result, class FmmTag>
class x2p_precompute;


template <class Result>
class x2p_precompute<Result, l2p_tag>
	: public fmm_operator<l2p_tag>
{
public:
	typedef Result result_t;

	template <class Operator>
	x2p_precompute(Operator const &op,
		cluster_tree<typename std::decay<Operator>::type::cluster_t> const &tree)
		: m_single_idx(tree.get_n_clusters())
	{
		auto const &list = tree.get_leaf_rec_indices();
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

private:
	std::vector<result_t> m_container;
	std::vector<size_t> m_single_idx;
};




template <class Result>
class x2p_precompute<Result, m2p_tag>
	: public fmm_operator<m2p_tag>
{
public:
	typedef interaction_lists::list_t list_t;
	typedef Result result_t;

	template <class Operator>
	x2p_precompute(Operator const &op,
		cluster_tree<typename std::decay<Operator>::type::cluster_t> const &tree,
		list_t const &list)
		: m_double_idx(tree.get_n_clusters(), tree.get_n_clusters())
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

	result_t const &operator()(size_t to, size_t from) const
	{
		return m_container[m_double_idx.coeff(to, from)];
	}

private:
	std::vector<result_t> m_container;
	Eigen::SparseMatrix<size_t> m_double_idx;
};



template <class Operator>
auto create_x2p_precompute(
	Operator const &op,
	cluster_tree<typename std::decay<Operator>::type::cluster_t> const &tree,
	interaction_lists::list_t const &list)
{
	typedef typename std::decay<Operator>::type operator_t;

	return x2p_precompute<
		typename operator_t::result_t,
		m2p_tag
	>(op, tree, list);
}

template <class Operator>
auto
create_x2p_precompute(
	Operator const &op, 
	cluster_tree<typename std::decay<Operator>::type::cluster_t> const &tree)
{
	typedef typename std::decay<Operator>::type operator_t;
	return x2p_precompute<
		typename operator_t::result_t,
		l2p_tag
	>(op, tree);
}


} // end of namespace fmm
} // namespace NiHu

#endif // LEAF_PRECOMPUTE_HPP_INCLUDED
