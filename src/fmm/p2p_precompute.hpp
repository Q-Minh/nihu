#ifndef P2P_PRECOMPUTE_HPP_INCLUDED
#define P2P_PRECOMPUTE_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "lists.hpp"
#include "util/matrix_traits.hpp"

#include <Eigen/SparseCore>

#include <vector>

namespace NiHu
{
namespace fmm
{

template <class Operator, class ClusterDerived>
Eigen::SparseMatrix<typename scalar<typename Operator::result_t>::type>
p2p_precompute(Operator const &op, cluster_tree<ClusterDerived> const &tree, interaction_lists::list_t const &list)
{
	typedef Operator operator_t;
	typedef typename operator_t::result_t result_t;
	typedef typename scalar<result_t>::type scalar_t;
	static size_t const rows = num_rows<result_t>::value;
	static size_t const cols = num_cols<result_t>::value;
	typedef Eigen::SparseMatrix<scalar_t> sparse_t;
	typedef Eigen::Triplet<scalar_t, size_t> triplet_t;

	std::vector<triplet_t> triplets;

	// determine number of entries in p2p sparse matrix
	size_t s = 0;
	for (unsigned to = 0; to < list.size(); ++to)
		for (auto from : list[to])
			s += tree[to].get_rec_node_idx().size() * rows *
			tree[from].get_src_node_idx().size() * cols;
	triplets.reserve(s);

	// compute entries and place them in p2p triplets
	for (unsigned to = 0; to < list.size(); ++to)
	{
		for (auto from : list[to])
		{
			for (auto i : tree[to].get_rec_node_idx())
			{
				for (auto j : tree[from].get_src_node_idx())
				{
					result_t mat = op(i, j);
					for (size_t ii = 0; ii < rows; ++ii)
						for (size_t jj = 0; jj < cols; ++jj)
							triplets.push_back(triplet_t(i*rows + ii, j*cols + jj, mat(ii, jj)));
				}
			}
		}
	}

	sparse_t mat(rows * tree.get_n_rec_nodes(), cols * tree.get_n_src_nodes());
	mat.setFromTriplets(triplets.begin(), triplets.end());
	return mat;
}

} // namespace fmm
} // namespace NiHu

#endif // P2P_PRECOMPUTE_HPP_INCLUDED
