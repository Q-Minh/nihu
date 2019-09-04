#ifndef FMM_X2X_PRECOMPUTE_HPP_INCLUDED
#define FMM_X2X_PRECOMPUTE_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "fmm_operator.hpp"
#include "lists.hpp"

#include "Eigen/SparseCore"

#include <vector>

namespace NiHu
{
namespace fmm
{

template <class Result, class ClusterDerived, class FmmTag>
class x2x_precompute
	: public fmm_operator<FmmTag>
{
public:
	typedef ClusterDerived cluster_t;
	typedef cluster_tree<cluster_t> cluster_tree_t;
	typedef interaction_lists::list_t list_t;
	typedef Result result_t;

	typedef std::vector<std::vector<result_t> > container_t;

	template <class Operator>
	x2x_precompute(Operator const &op, list_t const &list)
		: m_tree(op.get_tree())
		, m_indices(m_tree.get_n_clusters(), m_tree.get_n_clusters())
		, m_container(m_tree.get_n_levels())
	{
		typedef typename std::decay<Operator>::type operator_t;
		typedef Eigen::Triplet<size_t> triplet_t;
		std::vector<triplet_t> triplets;

		std::vector<std::vector<bool> > ready(m_tree.get_n_levels());

		for (size_t to = 0; to < list.size(); ++to)
		{
			for (auto from : list[to])
			{
				cluster_t const &cto = m_tree[to];
				cluster_t const &cfrom = m_tree[from];
				size_t level = cto.get_level();
				size_t idx = operator_t::operator_t::unique_idx(cto, cfrom);
				if (idx >= m_container[level].size())
				{
					m_container[level].resize(idx + 1);
					ready[level].resize(idx + 1, false);
				}
				if (!ready[level][idx])
				{
					m_container[level][idx] = op(to, from);
					ready[level][idx] = true;
				}
				triplets.push_back(triplet_t(int(to), int(from), idx));
			}
		}

		this->m_indices.setFromTriplets(triplets.begin(), triplets.end());
	}

	result_t const &operator()(size_t to, size_t from) const
	{
		size_t level = m_tree[to].get_level();
		size_t idx = m_indices.coeff(to, from);
		return m_container[level][idx];
	}

private:
	cluster_tree_t const &m_tree;
	Eigen::SparseMatrix<size_t> m_indices;
	container_t m_container;
};


template <class Operator>
auto
create_x2x_precompute(Operator const &op, typename interaction_lists::list_t const &list)
{
	typedef typename std::decay<Operator>::type operator_t;
	return x2x_precompute<
		typename operator_t::result_t,
		typename operator_t::cluster_t,
		typename operator_t::fmm_tag
	>(op, list);
}

} // namespace fmm
} // namespace NiHu

#endif // FMM_X2X_PRECOMPUTE_HPP_INCLUDED
