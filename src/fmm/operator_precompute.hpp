#ifndef OPERATOR_PRECOMPUTE_HPP_INCLUDED
#define OPERATOR_PRECOMPUTE_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "lists.hpp"

#include "Eigen/SparseCore"

#include <vector>

namespace NiHu
{
namespace fmm
{

template <class Operator>
class precompute
{
public:
	typedef Operator operator_t;
	typedef typename operator_t::cluster_t cluster_t;
	typedef cluster_tree<cluster_t> cluster_tree_t;
	typedef interaction_lists::list_t list_t;
	typedef typename operator_t::result_t result_t;

	typedef std::vector<std::vector<result_t> > container_t;

	precompute(operator_t const &op, cluster_tree_t const &tree, list_t const &list)
		: m_op(op)
		, m_tree(tree)
		, m_indices(tree.get_n_clusters(), tree.get_n_clusters())
		, m_container(tree.get_n_levels())
	{
		typedef Eigen::Triplet<size_t> triplet_t;
		std::vector<triplet_t> triplets;

		std::vector<std::vector<bool> > ready(tree.get_n_levels());

		for (size_t to = 0; to < list.size(); ++to)
		{
			for (auto from : list[to])
			{
				cluster_t const &cto = this->m_tree[to];
				cluster_t const &cfrom = this->m_tree[from];
				size_t level = cto.get_level();
				size_t idx = operator_t::unique_idx(cto, cfrom);
				if (idx >= this->m_container[level].size())
				{
					m_container[level].resize(idx + 1);
					ready[level].resize(idx + 1, false);
				}
				if (!ready[level][idx])
				{
					m_container[level][idx] = this->m_op(cto, cfrom);
					ready[level][idx] = true;
				}
				triplets.push_back(triplet_t(int(to), int(from), idx));
			}

		}

		this->m_indices.setFromTriplets(triplets.begin(), triplets.end());
	}

	result_t const &operator()(size_t to, size_t from) const
	{
		size_t level = this->m_tree[to].get_level();
		size_t idx = this->m_indices.coeff(to, from);
		return this->m_container[level][idx];
	}

private:
	operator_t const &m_op;
	cluster_tree_t const &m_tree;
	Eigen::SparseMatrix<size_t> m_indices;
	container_t m_container;
};

} // namespace fmm
} // namespace NiHu

#endif // OPERATOR_PRECOMPUTE_HPP_INCLUDED
