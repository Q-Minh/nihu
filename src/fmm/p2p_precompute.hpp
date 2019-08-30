#ifndef P2P_PRECOMPUTE_HPP_INCLUDED
#define P2P_PRECOMPUTE_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "fmm_operator.hpp"
#include "lists.hpp"
#include "util/eigen_utils.hpp"
#include "util/matrix_traits.hpp"

#include <Eigen/SparseCore>

#include <vector>
#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator, class ClusterDerived, bool isResultEigen = is_eigen<
	typename std::decay<Operator>::type::result_t>::value>
class p2p_precompute
	: public fmm_operator<p2p_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::result_t result_t;
	typedef typename scalar<result_t>::type scalar_t;
	static size_t const rows = num_rows<result_t>::value;
	static size_t const cols = num_cols<result_t>::value;
	typedef Eigen::SparseMatrix<scalar_t> sparse_t;
	typedef Eigen::Triplet<scalar_t, size_t> triplet_t;
	typedef ClusterDerived cluster_t;
	typedef cluster_tree<cluster_t> tree_t;

	p2p_precompute(Operator &&op, tree_t const &tree, interaction_lists::list_t const &list)
		: m_op(std::forward<Operator>(op))
		, m_tree(tree)
		, m_list(list)
	{
		compute_sparse_matrix();
	}

	void compute_sparse_matrix()
	{
		std::vector<triplet_t> triplets;

		// determine number of entries in p2p sparse matrix
		size_t s = 0;
		for (unsigned to = 0; to < m_list.size(); ++to)
			for (auto from : m_list[to])
				s += m_tree[to].get_rec_node_idx().size() * rows *
				m_tree[from].get_src_node_idx().size() * cols;
		triplets.reserve(s);

		// compute entries and place them in p2p triplets
		for (unsigned to = 0; to < m_list.size(); ++to)	// loop over receiver clusters
		{
			for (auto from : m_list[to])	// loop over source clusters
			{
				for (auto i : m_tree[to].get_rec_node_idx()) // loop over receiver nodes
				{
					for (auto j : m_tree[from].get_src_node_idx())	// loop over source nodes
					{
						result_t mat = m_op(i, j);
						for (size_t ii = 0; ii < rows; ++ii)	// loop over matrix rows
							for (size_t jj = 0; jj < cols; ++jj)	// loop over matrix cols
								triplets.push_back(triplet_t(i * rows + ii, j * cols + jj, mat(ii, jj)));
					}
				}
			}
		}

		m_mat.resize(rows * m_tree.get_n_rec_nodes(), cols * m_tree.get_n_src_nodes());
		m_mat.setFromTriplets(triplets.begin(), triplets.end());
	}

	sparse_t const &get_sparse_matrix() const
	{
		return m_mat;
	}

private:
	Operator m_op;
	tree_t const &m_tree;
	interaction_lists::list_t const &m_list;
	sparse_t m_mat;
};


template <class Operator, class ClusterDerived>
	class p2p_precompute<Operator, ClusterDerived, false>
	: public fmm_operator<p2p_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::result_t result_t;
	typedef typename scalar<result_t>::type scalar_t;
	static size_t const rows = num_rows<result_t>::value;
	static size_t const cols = num_cols<result_t>::value;
	typedef Eigen::SparseMatrix<scalar_t> sparse_t;
	typedef Eigen::Triplet<scalar_t, size_t> triplet_t;
	typedef ClusterDerived cluster_t;
	typedef cluster_tree<cluster_t> tree_t;

	p2p_precompute(Operator &&op, tree_t const &tree, interaction_lists::list_t const &list)
		: m_op(std::forward<Operator>(op))
		, m_tree(tree)
		, m_list(list)
	{
		compute_sparse_matrix();
	}

	void compute_sparse_matrix()
	{
		std::vector<triplet_t> triplets;

		// determine number of entries in p2p sparse matrix
		size_t s = 0;
		for (unsigned to = 0; to < m_list.size(); ++to)
			for (auto from : m_list[to])
				s += m_tree[to].get_rec_node_idx().size() * rows *
				m_tree[from].get_src_node_idx().size() * cols;
		triplets.reserve(s);

		// compute entries and place them in p2p triplets
		for (unsigned to = 0; to < m_list.size(); ++to)	// loop over receiver clusters
		{
			for (auto from : m_list[to])	// loop over source clusters
				for (auto i : m_tree[to].get_rec_node_idx()) // loop over receiver nodes
					for (auto j : m_tree[from].get_src_node_idx())	// loop over source nodes
						triplets.push_back(triplet_t(i, j, m_op(i, j)));
		}

		m_mat.resize(m_tree.get_n_rec_nodes(), m_tree.get_n_src_nodes());
		m_mat.setFromTriplets(triplets.begin(), triplets.end());
	}

	sparse_t const &get_sparse_matrix() const
	{
		return m_mat;
	}

private:
	Operator m_op;
	tree_t const &m_tree;
	interaction_lists::list_t const &m_list;
	sparse_t m_mat;
};



template <class Operator, class ClusterDerived>
auto create_p2p_precompute(Operator &&op,
	cluster_tree<ClusterDerived> const &tree,
	interaction_lists::list_t const &list)
{
	return p2p_precompute<Operator, ClusterDerived>(std::forward<Operator>(op), tree, list);
}


} // namespace fmm
} // namespace NiHu

#endif // P2P_PRECOMPUTE_HPP_INCLUDED
