/** \file lists.hpp
 * \brief declaration of class fmm::interaction_lists 
 */

#ifndef LISTS_HPP_INCLUDED
#define LISTS_HPP_INCLUDED

#include "cluster_tree.hpp"

#include <vector>
#include <iostream>
#include <iterator>

namespace NiHu
{
namespace fmm
{

/// \brief class storing the different interaction lists of a tree
class interaction_lists
{
public:
	enum {
		P2P = 0,
		P2M = 1,
		L2P = 2,
		P2L = 3,
		M2P = 4,
		M2M = 5,
		L2L = 6,
		M2L = 7,
		Near = 8
	};

	typedef std::vector<std::vector<size_t> > list_t;

private:
	static size_t const num_lists = 9;
	list_t lists[num_lists];
	static const bool p2l_allowed = false;

public:
	/** \brief constructor from cluster tree
	 * \tparam ClusterDerived the cluster type
	 * \param [in] tree the cluster tree
	 */
	template <class ClusterDerived>
	interaction_lists(cluster_tree<ClusterDerived> const &tree)
	{
		for (size_t i = 0; i < num_lists; ++i)
			lists[i].resize(tree.get_n_clusters());
		this->fill_lists(0, tree);
	}

	/** \brief return a selected list
	 * \param [in] idx the list index, interaction_lists::P2P for example
	 * \return the selected list
	 */
	list_t const &get_list(size_t idx) const
	{
		return this->lists[idx];
	}

private:
	/** \brief handle near field clusters of a leaf cluster
	 * \tparam ClusterDerived the cluster type
	 * \param [in] to index of the destination cluster
	 * \param [in] from index of the source cluster
	 * \param [in] tree the cluster tree
	 * This function handles the case when to is a leaf cluster and from is in
	 * its near field.
	 */
	template <class ClusterDerived>
	void add_adjacent_leaves(size_t to, size_t from, cluster_tree<ClusterDerived> const &tree)
	{
		if (tree[from].is_leaf())
		{
			if (tree[to].is_receiver() && tree[from].is_source())
				this->lists[P2P][to].push_back(from);
			if (tree[from].get_level() > tree[to].get_level())
				if (tree[from].is_receiver() && tree[to].is_source())
					this->lists[P2P][from].push_back(to);
		}
		else
		{
			for (auto pc : tree[from].get_children())
			{
				if (!interaction_lists::p2l_allowed 
					|| tree[to].get_bounding_box().is_adjacent(tree[pc].get_bounding_box()))
					add_adjacent_leaves(to, pc, tree);
				else
				{
					if (tree[to].is_receiver() && tree[pc].is_source())
						this->lists[M2P][to].push_back(pc);
					if (tree[pc].is_receiver() && tree[to].is_source())
						this->lists[P2L][pc].push_back(to);
				}
			}
		}
	}

	/** \brief fill the interaction lists of a cluster tree
	 * \tparam ClusterDerived the cluster type
	 * \param [in] root the index of the root cluster
	 * \param [in] tree the cluster tree
	 */
	template <class ClusterDerived>
	void fill_lists(size_t root, cluster_tree<ClusterDerived> const &tree)
	{
		// the root is in itself's near field
		if (root == 0)
			this->lists[Near][root].push_back(root);
		else
		{
			// find the parent's near field
			for (auto ppn : this->lists[Near][tree[root].get_parent()])
			{
				// traverse the parent's near field's children
				for (auto p : tree[ppn].get_children())
				{
					// adjacent cousin : near field
					if (tree[root].get_bounding_box().is_adjacent(tree[p].get_bounding_box()))
						this->lists[Near][root].push_back(p);
					else // non-adjacent cousin: interaction list
						if (tree[root].is_receiver() && tree[p].is_source())
							this->lists[M2L][root].push_back(p);
				}
			}
		}

		if (tree[root].is_leaf())
			for (auto pc : this->lists[Near][root])
				add_adjacent_leaves(root, pc, tree);
		else
		{
			for (auto pc : tree[root].get_children())
			{
				if (tree[root].get_level() >= 2)
				{
					if (tree[root].is_source())
						this->lists[M2M][root].push_back(pc);
					if (tree[root].is_receiver())
						this->lists[L2L][pc].push_back(root);
				}

				fill_lists(pc, tree);
			}
		}
	}

public:
	/** \brief print debug information of the lists to a stream
	 * \param [in] os the output stream
	 */
	void print_debug(std::ostream &os = std::cout)
	{
		std::ostream_iterator<size_t> out_it(os, ", ");

		for (size_t c = 0; c < this->lists[0].size(); ++c)
		{
			std::cout << c << std::endl;
			os << "(M2L): ";
			std::copy(this->lists[M2L][c].begin(), this->lists[M2L][c].end(), out_it);
			os << std::endl;

			os << "(M2P): ";
			std::copy(this->lists[M2P][c].begin(), this->lists[M2P][c].end(), out_it);
			os << std::endl;

			os << "(P2L): ";
			std::copy(this->lists[P2L][c].begin(), this->lists[P2L][c].end(), out_it);
			os << std::endl;

			os << "(P2P): ";
			std::copy(this->lists[P2P][c].begin(), this->lists[P2P][c].end(), out_it);
			os << std::endl;
		}
	}
};

} // end of namespace fmm
} // namespace NiHu

#endif // LISTS_HPP_INCLUDED
