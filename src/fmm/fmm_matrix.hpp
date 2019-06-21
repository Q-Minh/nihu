/// \file fmm_matrix.hpp
/// \brief definition of class fmm::fmm_matrix

#ifndef FMM_MATRIX_HPP_INCLUDED
#define FMM_MATRIX_HPP_INCLUDED

#include <omp.h>

#include <vector>

#include "cluster_tree.hpp"
#include "fmm_timer.hpp"
#include "lists.hpp"
#include "util/matrix_traits.hpp"

namespace NiHu
{
namespace fmm
{

/// \brief matrix representation of the FMM method
/// \tparam P2P the P2P operator's type (sparse matrix)
/// \tparam P2M the P2M operator's type
/// \tparam P2L the P2L operator's type
/// \tparam M2P the M2P operator's type
/// \tparam L2P the L2P operator's type
/// \tparam M2M the M2M operator's type
/// \tparam M2L the M2L operator's type
/// \tparam L2L the L2L operator's type
template <
	class P2P,
	class P2M, class P2L, class M2P, class L2P,
	class M2M, class L2L, class M2L
>
class fmm_matrix
{
public:
	typedef typename std::decay<P2P>::type p2p_t;
	typedef typename std::decay<P2M>::type p2m_t;
	typedef typename std::decay<P2L>::type p2l_t;
	typedef typename std::decay<M2P>::type m2p_t;
	typedef typename std::decay<L2P>::type l2p_t;
	typedef typename std::decay<M2M>::type m2m_t;
	typedef typename std::decay<L2L>::type l2l_t;
	typedef typename std::decay<M2L>::type m2l_t;

	typedef typename m2l_t::cluster_t cluster_t;
	typedef typename cluster_t::multipole_t multipole_t;
	typedef typename cluster_t::local_t local_t;

	typedef typename scalar<p2p_t>::type scalar_t;
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> excitation_t;
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> response_t;

	typedef cluster_tree<cluster_t> cluster_tree_t;

	/// \brief number of DOF for a source node in the mesh 
	static size_t const num_dof_per_src = p2m_t::num_dof_per_src;
	/// \brief number of DOF for a receiver node in the mesh 
	static size_t const num_dof_per_rec = l2p_t::num_dof_per_rec;

	/// \brief constructor from operator instances
	/// \param [in] p2p the P2P operator
	/// \param [in] p2m the P2M operator
	/// \param [in] p2l the P2L operator
	/// \param [in] m2p the M2P operator
	/// \param [in] l2p the L2P operator
	/// \param [in] m2m the M2M operator
	/// \param [in] l2l the L2L operator
	/// \param [in] m2l the M2L operator
	/// \param [in] tree the cluster tree
	/// \param [in] lists the interaction lists
	fmm_matrix(P2P &&p2p,
		P2M &&p2m, P2L &&p2l, M2P &&m2p, L2P &&l2p,
		M2M &&m2m, L2L &&l2l, M2L &&m2l,
		cluster_tree_t const &tree,
		interaction_lists const &lists)
		: m_p2p(std::forward<P2P>(p2p))
		, m_p2m(std::forward<P2M>(p2m))
		, m_p2l(std::forward<P2L>(p2l))
		, m_m2p(std::forward<M2P>(m2p))
		, m_l2p(std::forward<L2P>(l2p))
		, m_m2m(std::forward<M2M>(m2m))
		, m_l2l(std::forward<L2L>(l2l))
		, m_m2l(std::forward<M2L>(m2l))
		, m_tree(tree)
		, m_lists(lists)
		, m_timer(tree.get_n_levels())
		, m_rhs_segments(tree.get_n_clusters())
		, m_lhs_segments(tree.get_n_clusters())
	{
		for (auto c : m_tree.get_leaf_src_indices())
		{
			cluster_t const &clus = m_tree[c];
			m_rhs_segments[c].resize(clus.get_n_src_nodes() * num_dof_per_src);
		}
		for (auto c : m_tree.get_leaf_rec_indices())
		{
			cluster_t const &clus = m_tree[c];
			m_lhs_segments[c].resize(clus.get_n_rec_nodes() * num_dof_per_rec);
		}
	}

	/// \brief return number of rows of the matrix
	/// \return number of rows
	size_t rows() const
	{
		return m_p2p.rows();
	}

	/// \brief return number of columns of the matrix
	/// \return number of columns
	size_t cols() const
	{
		return m_p2p.cols();
	}

	/// \brief cluster-continuous reordering of the excitation data
	/// \details after reordering, the data associated with a specific cluster can be
	/// reached as a continuous block (segment) of the internally stored excitation vector
	/// \param [in] rhs the right hand side vector
	template <class Derived>
	void reorder_excitation(Eigen::MatrixBase<Derived> const &rhs)
	{
		for (auto c : m_tree.get_leaf_src_indices())
		{
			cluster_t const &clus = m_tree[c];
			for (size_t i = 0; i < clus.get_n_src_nodes(); ++i)
			{
				size_t ii = clus.get_src_node_idx()[i];
				m_rhs_segments[c].segment(i*num_dof_per_src, num_dof_per_src) =
					rhs.segment(ii*num_dof_per_src, num_dof_per_src);
			}
		}
	}

	void upward_pass_dfs_rec(std::vector<multipole_t> &multipoles, size_t root)
	{
		for (auto c : m_tree[root].get_children())
		{
			if (m_tree[c].is_source())
			{
				upward_pass_dfs_rec(multipoles, c);
				multipoles[root] += m_m2m(root, c) * multipoles[c];
			}
		}
	}

	void downward_pass_dfs_rec(std::vector<local_t> &locals,
		std::vector<multipole_t> const &multipoles, size_t to)
	{
		if (!m_tree[to].is_receiver())
			return;

		for (auto from : m_lists.get_list(interaction_lists::M2L)[to])
			locals[to] += m_m2l(to, from) * multipoles[from];

		for (auto from : m_lists.get_list(interaction_lists::L2L)[to])
			locals[to] += m_l2l(to, from) * locals[from];

		for (auto c : m_tree[to].get_children())
			downward_pass_dfs_rec(locals, multipoles, c);
	}

	void upward_pass_dfs(std::vector<multipole_t> &multipoles)
	{
		size_t a = m_tree.level_begin(2);
		size_t b = m_tree.level_end(2);
#ifdef PARALLEL
#pragma omp parallel for
#endif
		for (int to = a; to < b; ++to)
			upward_pass_dfs_rec(multipoles, to);
	}

	void downward_pass_dfs(std::vector<local_t> &locals,
		std::vector<multipole_t> const &multipoles)
	{
		size_t a = m_tree.level_begin(2);
		size_t b = m_tree.level_end(2);
#ifdef PARALLEL
#pragma omp parallel for
#endif
		for (int to = int(a); to < int(b); ++to)
			downward_pass_dfs_rec(locals, multipoles, to);
	}

	void upward_pass_bfs(std::vector<multipole_t> &multipoles)
	{
		size_t nLevels = m_tree.get_n_levels();

		// compute upward pass
		for (size_t iLevel = nLevels - 2; iLevel >= 2; --iLevel)
		{
			m_timer.tic();
#ifdef PARALLEL
#pragma omp parallel for
#endif
			for (int to = int(m_tree.level_begin(iLevel)); to < int(m_tree.level_end(iLevel)); ++to)
			{
				for (auto from : m_tree[to].get_children())
				{
					if (!m_tree[from].is_source())
						continue;
					multipoles[to] += m_m2m(to, from) * multipoles[from];
				}
			}
#pragma omp barrier
			m_timer.toc(iLevel, fmm_timer::M2M);
		}
	}

	void downward_pass_bfs(std::vector<local_t> &locals,
		std::vector<multipole_t> const &multipoles)
	{
		size_t nLevels = m_tree.get_n_levels();

		for (unsigned iLevel = 2; iLevel < nLevels; ++iLevel)
		{
			size_t a = m_tree.level_begin(iLevel);
			size_t b = m_tree.level_end(iLevel);

			m_timer.tic();
#ifdef PARALLEL
#pragma omp parallel for
#endif
			for (int to = int(a); to < int(b); ++to)
				for (size_t from : m_lists.get_list(interaction_lists::M2L)[to])
					locals[to] += m_m2l(to, from) * multipoles[from];
#pragma omp barrier
			m_timer.toc(iLevel, fmm_timer::M2L);

			// no L2L needed at highest level
			if (iLevel == 2)
				continue;

			m_timer.tic();
#ifdef PARALLEL
#pragma omp parallel for
#endif
			for (int to = int(a); to < int(b); ++to)
			{
				if (!m_tree[to].is_receiver())
					continue;
				size_t from = m_tree[to].get_parent();
				locals[to] += m_l2l(to, from) * locals[from];
			}
#pragma omp barrier
			m_timer.toc(iLevel, fmm_timer::L2L);

		}
	}


	/// \brief cluster-continuous inverse-reordering of the response data
	/// \details after reordering, the response data contains the response in the same
	/// order as defined by the excitation
	/// \param [out] lhs the response vector
	template <class Derived>
	void reorder_response(Eigen::MatrixBase<Derived> &lhs) const
	{
		for (auto c : m_tree.get_leaf_rec_indices())
		{
			cluster_t const &clus = m_tree[c];
			for (size_t i = 0; i < clus.get_n_rec_nodes(); ++i)
			{
				size_t ii = clus.get_rec_node_idx()[i];
				lhs.segment(ii*num_dof_per_rec, num_dof_per_rec)
					+= m_lhs_segments[c].segment(i*num_dof_per_rec, num_dof_per_rec);
			}
		}
	}

	/// \brief matrix vector multiplication
	/// \param [in] rhs the excitation vector
	/// \return the result of matrix_vector multiplication 
	template <class ExcType>
	response_t operator *(ExcType const &rhs)
	{
		// compute P2P interactions
		m_timer.tic();
		response_t lhs = m_p2p * rhs;
		m_timer.toc(0, fmm_timer::P2P);

		// instantiate locals and multipoles
		size_t n_clusters = m_tree.get_n_clusters();
		std::vector<multipole_t> multipoles(n_clusters);
		std::vector<local_t> locals(n_clusters);
		for (size_t c = 0; c < n_clusters; ++c)
		{
			if (m_tree[c].get_level() < 2)
				continue;
			if (m_tree[c].is_source())
				multipoles[c] = m_tree[c].zero_multipole();
			if (m_tree[c].is_receiver())
				locals[c] = m_tree[c].zero_local();
		}

		// read reordered excitation into cluster data
		reorder_excitation(rhs);

		// compute P2L interactions
		m_timer.tic();
		auto const &P2Llist = m_lists.get_list(interaction_lists::P2L);
		for (auto to : m_tree.get_leaf_rec_indices())
			for (auto from : P2Llist[to])
				locals[to] += m_p2l(to, from) * m_rhs_segments[from];
		m_timer.toc(0, fmm_timer::P2L);

		// compute P2M interactions
		m_timer.tic();
		auto const &src_idx = m_tree.get_leaf_src_indices();
#ifdef PARALLEL
#pragma omp parallel for
#endif
		for (int i = 0; i < int(src_idx.size()); ++i)
		{
			size_t to = src_idx[i];
			multipoles[to] += m_p2m(to) * m_rhs_segments[to];
		}
#pragma omp barrier
		m_timer.toc(0, fmm_timer::P2M);

#if defined BFS
		upward_pass_bfs(multipoles);
		downward_pass_bfs(locals, multipoles);
#elif defined DFS
		upward_pass_dfs(multipoles);
		downward_pass_dfs(locals, multipoles);
#else
#	error You need to explicitly define BFS or DFS tree traversing algorithm
#endif

		// compute L2P interactions
		m_timer.tic();
		auto const &rec_idx = m_tree.get_leaf_rec_indices();
#ifdef PARALLEL
#pragma omp parallel for
#endif
		for (int i = 0; i < int(rec_idx.size()); ++i)
		{
			size_t to = rec_idx[i];
			m_lhs_segments[to] = m_l2p(to) * locals[to];
		}
#pragma omp barrier
		m_timer.toc(0, fmm_timer::L2P);

		// compute M2P interactions
		m_timer.tic();
		for (auto to : m_tree.get_leaf_rec_indices())
			for (auto from : m_lists.get_list(interaction_lists::M2P)[to])
				m_lhs_segments[to] += m_m2p(to, from) * multipoles[from];
		m_timer.toc(0, fmm_timer::M2P);

		// distribute reordered response
		reorder_response(lhs);

		return lhs;
	}

	fmm_timer const &get_timer() const
	{
		return m_timer;
	}

	Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> get_diagonal() const
	{
		return m_p2p.diagonal();
	}

private:
	P2P m_p2p;
	P2M m_p2m;
	P2L m_p2l;
	M2P m_m2p;
	L2P m_l2p;
	M2M m_m2m;
	L2L m_l2l;
	M2L m_m2l;
	cluster_tree_t const &m_tree;
	interaction_lists const &m_lists;
	fmm_timer m_timer;
	std::vector<excitation_t> m_rhs_segments;
	std::vector<response_t> m_lhs_segments;
};


/// \brief factory function to create an fmm_matrix object
/// \param [in] p2p the P2P operator object
/// \param [in] p2m the P2M operator object
/// \param [in] p2l the P2L operator object
/// \param [in] m2p the M2P operator object
/// \param [in] l2p the L2P operator object
/// \param [in] m2m the M2M operator object
/// \param [in] l2l the L2L operator object
/// \param [in] m2l the M2L operator object
/// \param [in] tree the cluster tree
/// \param [in] lists the interaction lists
template <class P2P, class P2M, class P2L, class M2P, class L2P,
	class M2M, class L2L, class M2L, class Cluster>
	fmm_matrix<P2P, P2M, P2L, M2P, L2P, M2M, L2L, M2L>
	create_fmm_matrix(
		P2P &&p2p,
		P2M &&p2m, P2L &&p2l, M2P &&m2p, L2P &&l2p,
		M2M &&m2m, L2L &&l2l, M2L &&m2l,
		cluster_tree<Cluster> const &tree,
		interaction_lists const &lists
	)
{
	return fmm_matrix<P2P, P2M, P2L, M2P, L2P, M2M, L2L, M2L>(
		std::forward<P2P>(p2p),
		std::forward<P2M>(p2m),
		std::forward<P2L>(p2l),
		std::forward<M2P>(m2p),
		std::forward<L2P>(l2p),
		std::forward<M2M>(m2m),
		std::forward<L2L>(l2l),
		std::forward<M2L>(m2l),
		tree,
		lists);
}

} // end of namespace fmm
} // namespace NiHu

#endif
