/** 
 * @file cluster_tree.hpp
 * @brief Implementation of class @ref NiHu::fmm::cluster_tree
 * @ingroup fmm_clus
 */

#ifndef FMM_CLUSTER_TREE_HPP_INCLUDED
#define FMM_CLUSTER_TREE_HPP_INCLUDED

#include "cluster.hpp"
#include "divide.hpp"
#include "util/misc.hpp"

#include "Eigen/Dense"
#include "Eigen/StdVector"

#include <algorithm> 	// std::min_element
#include <cstddef>
#include <iostream>
#include <numeric> 		// std::iota
#include <stdexcept>	// std::logic_error

namespace NiHu
{
namespace fmm
{

/** 
 * @brief Class representing a cluster tree
 * @tparam Dim Space dimension
 */
template <class ClusterDerived>
class cluster_tree
{
public:
	typedef ClusterDerived cluster_t;

	static size_t const dimension = cluster_t::dimension;
	typedef typename cluster_t::bounding_box_t bounding_box_t;
	typedef typename cluster_t::location_t location_t;
	typedef typename cluster_t::idx_list_t idx_list_t;

	/** @brief Type of the cluster vector */
	typedef std::vector<
		cluster_t, Eigen::aligned_allocator<cluster_t>
	> cluster_vector_t;

	typedef typename cluster_vector_t::const_iterator iterator_t;

public:
	/**
	 * @brief Create a cluster tree 
	 * @tparam It Iterator type
	 * @tparam DivideDerived Cluster division method 
	 * @param[in] src_begin Begin iterator to sources / receivers
	 * @param[in] src_end End iterator to sources / receivers
	 */
	template <class It, class DivideDerived>
	cluster_tree(It src_begin, It src_end, divide_base<DivideDerived> const& divide)
		: cluster_tree(src_begin, src_end, src_begin, src_end, divide)
	{
	}

	/**
	 * @brief Create a cluster tree 
	 * @tparam It1 Source iterator type
	 * @tparam It2 Receiver iterator type
	 * @tparam DivideDerived Cluster division method 
	 * @param[in] src_begin Begin iterator to sources
	 * @param[in] src_end End iterator to sources 
	 * @param[in] rec_begin Begin iterator to receivers 
	 * @param[in] rec_end End iterator to receivers
	 */
	template <class It1, class It2, class DivideDerived>
	cluster_tree(It1 src_begin, It1 src_end, It2 rec_begin, It2 rec_end, divide_base<DivideDerived> const& divide)
		: n_src(src_end - src_begin)
		, n_rec(rec_end - rec_begin)
	{
		if (this->n_src + this->n_rec == 0)
			throw std::logic_error("Cannot build cluster_tree from zero nodes");

		cluster_t root_cluster;
		root_cluster.set_level(0);

		// compute bounding box of source and receiver nodes
		typedef Eigen::Matrix<double, dimension, Eigen::Dynamic> nodes_t;
		nodes_t nodes(dimension, n_src + n_rec);
		for (size_t i = 0; i < n_src; ++i)
			nodes.col(i) = src_begin[i];
		for (size_t i = 0; i < n_rec; ++i)
			nodes.col(n_src + i) = rec_begin[i];
		bounding_box_t bb(nodes);
		root_cluster.set_bounding_box(bb);

		// compute node index vectors of root cluster
		idx_list_t src_node_idx(n_src);
		std::iota(src_node_idx.begin(), src_node_idx.end(), 0);
		root_cluster.set_src_node_idx(src_node_idx);

		idx_list_t rec_node_idx(n_rec);
		std::iota(rec_node_idx.begin(), rec_node_idx.end(), 0);
		root_cluster.set_rec_node_idx(rec_node_idx);

		this->clusters.push_back(root_cluster);

		// BFS traverse
		for (size_t idx = 0; idx < clusters.size(); ++idx)
		{
			if (divide(clusters[idx]))
			{
				// get parent bounding box
				auto par_bb = clusters[idx].get_bounding_box();

				// copy of parent source and receiver node indices
				auto par_src = clusters[idx].get_src_node_idx();
				auto par_rec = clusters[idx].get_rec_node_idx();

				// traverse possible children
				for (size_t i = 0; i < (1 << dimension); ++i)
				{
					auto child_bb = par_bb.get_child(i);

					// move source indices to child if contained
					idx_list_t child_src(par_src.size());
					auto q_src = move_if(par_src.begin(), par_src.end(),
						child_src.begin(), [&](size_t u) {
						return bounding_box_t::dist2idx(src_begin[u], par_bb.get_center()) == i;
					});
					child_src.resize(std::distance(q_src, par_src.end()));
					par_src.erase(q_src, par_src.end());

					// move receiver indices to child if contained
					idx_list_t child_rec(par_rec.size());
					auto q_rec = move_if(par_rec.begin(), par_rec.end(),
						child_rec.begin(), [&](size_t u) {
						return bounding_box_t::dist2idx(rec_begin[u], par_bb.get_center()) == i;
					});
					child_rec.resize(std::distance(q_rec, par_rec.end()));
					par_rec.erase(q_rec, par_rec.end());

					// if child is empty, don't add to tree
					if (child_src.empty() && child_rec.empty())
						continue;

					cluster_t child;
					child.set_parent(idx);
					child.set_bounding_box(child_bb);
					child.set_src_node_idx(child_src);
					child.set_rec_node_idx(child_rec);
					child.set_level(clusters[idx].get_level() + 1);
					clusters[idx].add_child(clusters.size());
					clusters.push_back(child);
				}
				if (!par_src.empty() || !par_rec.empty())
					throw std::runtime_error("Remaining nodes in divided cluster");
			}
			else
			{
				this->leaf_indices.push_back(idx);
				if (clusters[idx].is_source())
					this->leaf_src_indices.push_back(idx);
				if (clusters[idx].is_receiver())
					this->leaf_rec_indices.push_back(idx);
			}
		}

		// build levels structure
		this->levels.push_back(0);
		for (size_t i = 0; i < this->clusters.size() - 1; ++i)
			if (this->clusters[i].get_level() != this->clusters[i + 1].get_level())
				this->levels.push_back(i + 1);
		this->levels.push_back(this->clusters.size());

	}

	size_t get_n_src_nodes() const
	{
		return this->n_src;
	}

	size_t get_n_rec_nodes() const
	{
		return this->n_rec;
	}

	/**\brief return vector of leaf cluster indices
	 * \return vector of leaf cluster indices
	 */
	std::vector<size_t> const &get_leaf_indices() const
	{
		return this->leaf_indices;
	}

	std::vector<size_t> const &get_leaf_src_indices() const
	{
		return this->leaf_src_indices;
	}

	std::vector<size_t> const &get_leaf_rec_indices() const
	{
		return this->leaf_rec_indices;
	}

	/**\brief return number of leaf clusters
	 * \return number of leaf clusters
	 */
	size_t get_n_leaves() const
	{
		return this->leaf_indices.size();
	}

	/**\brief assemble and return source leaf index vectors
	 * \return vector of vectors containing leaf source indices
	 */
	std::vector<idx_list_t> get_leaf_src_index_vectors() const
	{
		std::vector<idx_list_t> indices(get_n_leaves());
		auto v = get_leaf_indices();
		for (size_t i = 0; i < v.size(); ++i)
			indices[i] = this->clusters[v[i]].get_src_node_idx();
		return indices;
	}

	/**\brief index operator returning the idx-th cluster
	 * \param [in] idx the cluster index
	 * \return the idx-th cluster
	 */
	cluster_t const &operator[](size_t idx) const
	{
		return clusters[idx];
	}

	/**\brief index operator returning the idx-th cluster
	 * \param [in] idx the cluster index
	 * \return the idx-th cluster
	 */
	cluster_t &operator[](size_t idx)
	{
		return clusters[idx];
	}

	/**\brief return root bounding box diameter
	 * \return root bounding box diameter
	 */
	double get_root_diameter() const
	{
		return clusters[0].get_bounding_box().get_diameter();
	}

	/**\brief return number of levels
	 * \return number of levels
	 */
	size_t get_n_levels() const
	{
		return this->levels.size() - 1;
	}

	/**\brief return begin iterator to idx-th level clusters
	 * \param [in] idx the level index
	 * \return begin index of idx-th level clusters
	 */
	size_t level_begin(size_t idx) const
	{
		return this->levels[idx];
	}

	/**\brief return end iterator to idx-th level clusters
	 * \param [in] idx the level index
	 * \return end index of idx-th level clusters
	 */
	size_t level_end(size_t idx) const
	{
		return this->levels[idx + 1];
	}

	/**\brief return number of clusters
	 * \return number of clusters
	 */
	size_t get_n_clusters() const
	{
		return this->clusters.size();
	}

	iterator_t begin() const
	{
		return clusters.begin();
	}

	iterator_t end() const
	{
		return clusters.end();
	}

	/**\brief print debug information to output strean
	 * \param [in] os the output stream
	 */
	void print_debug(std::ostream &os = std::cout) const
	{
		os << "#Levels: " << this->get_n_levels() << std::endl;
		for (size_t l = 0; l < this->get_n_levels(); ++l)
			os << "Level #" << l << ": " << level_end(l) - level_begin(l) << " clusters" << std::endl;
		os << "#Source Nodes: " << this->get_n_src_nodes() << std::endl;
		os << "#Receiver Nodes: " << this->get_n_rec_nodes() << std::endl;
		os << "Root diameter: " << this->get_root_diameter() << std::endl;
	}

	idx_list_t get_src_ids() const
	{
		return this->src_node_ids;
	}

	idx_list_t const &get_rec_ids() const
	{
		return this->rec_node_id;
	}

private:
	size_t n_src, n_rec;

	cluster_vector_t clusters;

	std::vector<size_t> leaf_indices;
	std::vector<size_t> leaf_src_indices;
	std::vector<size_t> leaf_rec_indices;
	std::vector<size_t> levels;
};

template <class ClusterDerived>
size_t const cluster_tree<ClusterDerived>::dimension;


template <class C>
std::ostream &operator<<(std::ostream &os, cluster_tree<C> const &ct)
{
	ct.print_debug(os);
	return os;
}

} // end of namespace fmm
} // end of namespace NiHu

#endif /* FMM_CLUSTER_TREE_HPP_INCLUDED */
