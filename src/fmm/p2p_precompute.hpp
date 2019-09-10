/**
 * @file p2p_precompute.hpp 
 * @brief Pre-computation of P2P operators for acceleration 
 * @ingroup fmm_ops
 */

#ifndef P2P_PRECOMPUTE_HPP_INCLUDED
#define P2P_PRECOMPUTE_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "fmm_operator.hpp"
#include "lists.hpp"
#include "util/eigen_utils.hpp"
#include "util/matrix_traits.hpp"

#include <Eigen/SparseCore>

#include <chrono>
#include <omp.h>
#include <type_traits>
#include <utility>
#include <vector>

namespace NiHu
{
namespace fmm
{

namespace internal
{

template <class Operator, class ClusterDerived, bool isResultEigen = is_eigen<
	typename std::decay<Operator>::type::result_t>::value>
class sparse_computer;


template <class Operator, class ClusterDerived>
class sparse_computer<Operator, ClusterDerived, true>
{
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::result_t result_t;
	typedef typename scalar<result_t>::type scalar_t;
	static size_t const rows = num_rows<result_t>::value;
	static size_t const cols = num_cols<result_t>::value;
	typedef Eigen::Triplet<scalar_t, size_t> triplet_t;
	typedef Eigen::SparseMatrix<scalar_t> sparse_t;
	typedef ClusterDerived cluster_t;
	typedef cluster_tree<cluster_t> tree_t;

public:

	static sparse_t eval(Operator &&op, tree_t const &tree, interaction_lists::list_t const &list, size_t &assembly_time)
	{
		auto tstart = std::chrono::steady_clock::now();
	
		std::vector<triplet_t> triplets;
		
		std::vector<std::pair<size_t, size_t> > indices;

		// determine number of src receiver pairs
		size_t s = 0;
		for (unsigned to = 0; to < list.size(); ++to)
			for (auto from : list[to])
				s += tree[to].get_rec_node_idx().size() * 
				tree[from].get_src_node_idx().size();
				
		indices.reserve(s);
		
		for (unsigned to = 0; to < list.size(); ++to)
			for (auto from : list[to])
				for (auto i : tree[to].get_rec_node_idx()) // loop over receiver nodes
					for (auto j : tree[from].get_src_node_idx())	// loop over source nodes
						indices.push_back(std::make_pair(i, j));
						
		triplets.resize(s * rows * cols);
		
#ifdef PARALLEL
#pragma omp parallel for
#endif
		for (int isrcrec = 0; isrcrec < s; ++isrcrec)
		{
			size_t i = indices[isrcrec].first;
			size_t j = indices[isrcrec].second;
			result_t mat = op(i, j);
			for (size_t ii = 0; ii < rows; ++ii)	// loop over matrix rows
				for (size_t jj = 0; jj < cols; ++jj)	// loop over matrix cols
					triplets[isrcrec*rows*cols + ii*cols + jj] = triplet_t(i * rows + ii, j * cols + jj, mat(ii, jj));
		}
#ifdef PARALLEL
#pragma omp barrier
#endif

		sparse_t mat(rows * tree.get_n_rec_nodes(), cols * tree.get_n_src_nodes());
		mat.setFromTriplets(triplets.begin(), triplets.end());
		
		auto tend = std::chrono::steady_clock::now();
		assembly_time = std::chrono::duration_cast<std::chrono::microseconds>(tend - tstart).count();
		
		return mat;
	}
};

template <class Operator, class ClusterDerived>
class sparse_computer<Operator, ClusterDerived, false>
{
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::result_t result_t;
	typedef typename scalar<result_t>::type scalar_t;
	typedef Eigen::Triplet<scalar_t, size_t> triplet_t;
	typedef Eigen::SparseMatrix<scalar_t> sparse_t;
	typedef ClusterDerived cluster_t;
	typedef cluster_tree<cluster_t> tree_t;

public:

	static sparse_t eval(Operator &&op, tree_t const &tree, interaction_lists::list_t const &list, size_t &assembly_time)
	{
		auto tstart = std::chrono::steady_clock::now();
	
		std::vector<triplet_t> triplets;
		
		std::vector<std::pair<size_t, size_t> > indices;

		// determine number of src receiver pairs
		size_t s = 0;
		for (unsigned to = 0; to < list.size(); ++to)
			for (auto from : list[to])
				s += tree[to].get_rec_node_idx().size() * 
				tree[from].get_src_node_idx().size();
				
		indices.reserve(s);
		
		for (unsigned to = 0; to < list.size(); ++to)
			for (auto from : list[to])
				for (auto i : tree[to].get_rec_node_idx()) // loop over receiver nodes
					for (auto j : tree[from].get_src_node_idx())	// loop over source nodes
						indices.push_back(std::make_pair(i, j));
						
		triplets.resize(s);
		
#ifdef PARALLEL
#pragma omp parallel for
#endif
		for (int isrcrec = 0; isrcrec < s; ++isrcrec)
		{
			size_t i = indices[isrcrec].first;
			size_t j = indices[isrcrec].second;
			triplets[isrcrec] = triplet_t(i, j, op(i, j));
		}
#ifdef PARALLEL
#pragma omp barrier
#endif

		sparse_t mat(tree.get_n_rec_nodes(), tree.get_n_src_nodes());
		mat.setFromTriplets(triplets.begin(), triplets.end());
		
		auto tend = std::chrono::steady_clock::now();
		assembly_time = std::chrono::duration_cast<std::chrono::microseconds>(tend - tstart).count();
		
		return mat;
	}
};


} // end of namespace internal

template <class Scalar, size_t NumDofPerSrc, size_t NumDofPerRec>
class p2p_precompute
	: public fmm_operator<p2p_tag>
{
public:
	typedef Scalar scalar_t;
	typedef Eigen::SparseMatrix<scalar_t> sparse_t;

	static size_t const num_dof_per_src = NumDofPerSrc;
	static size_t const num_dof_per_rec = NumDofPerRec;

	template <class Operator, class ClusterDerived>
	p2p_precompute(Operator &&op, cluster_tree<ClusterDerived> const &tree, interaction_lists::list_t const &list)
		: m_mat(internal::sparse_computer<Operator, ClusterDerived>::eval(std::forward<Operator>(op), tree, list, m_assembly_time))
	{
	}

	sparse_t const &get_sparse_matrix() const
	{
		return m_mat;
	}
	
	size_t get_assembly_time() const
	{
		return m_assembly_time;
	}

private:
	sparse_t m_mat;
	size_t m_assembly_time;
};


template <class Operator, class ClusterDerived>
auto create_p2p_precompute(Operator &&op,
	cluster_tree<ClusterDerived> const &tree,
	interaction_lists::list_t const &list)
{
	typedef typename std::decay<Operator>::type operator_t;
	return p2p_precompute<
		typename scalar<typename operator_t::result_t>::type,
		operator_t::num_dof_per_src,
		operator_t::num_dof_per_rec
	>(std::forward<Operator>(op), tree, list);
}


} // end of namespace fmm
} // end of namespace NiHu

#endif /* P2P_PRECOMPUTE_HPP_INCLUDED */
