/** \file cluster.hpp
 * \brief implementation of class NiHu::fmm::cluster_base
 */
#ifndef CLUSTER_HPP_INCLUDED
#define CLUSTER_HPP_INCLUDED

#include "bounding_box.hpp"

#include <algorithm>	// std::find, std::copy
#include <iostream>		// std::ostream
#include <iterator>		// std::ostream_iterator
#include <vector>		// idx_list
#include <stdexcept>	// std::logic_error

namespace NiHu
{
namespace fmm
{
/** \brief CRTP traits structure of a cluster
 * \tparam Derived the CRTP derived class
 *
 * The traits structure needs to define
 * local_t,
 * multipole_t,
 * dimension
 */
template <class Derived>
struct cluster_traits;

/** \brief CRTP base class of clusters
 * \tparam Derived CRTP derived class
 */
template <class Derived>
class cluster_base
{
public:
	/** \brief the traits structure */
	typedef cluster_traits<Derived> traits_t;

	/** \brief the space dimension */
	static size_t const dimension = traits_t::dimension;
	/** \brief the multipole type */
	typedef typename traits_t::multipole_t multipole_t;
	/** \brief the local type */
	typedef typename traits_t::local_t local_t;

	/** \brief the bounding box type */
	typedef bounding_box<dimension, double> bounding_box_t;
	/** \brief the location type */
	typedef typename bounding_box_t::location_t location_t;
	/** \brief index list type for storing children and contained nodes */
	typedef std::vector<size_t> idx_list_t;

public:
	/** \brief crtp helper function */
	Derived const &derived() const
	{
		return *(static_cast<Derived const *>(this));
	}

	/** \brief crtp helper function */
	Derived &derived()
	{
		return *(static_cast<Derived *>(this));
	}

	/** \brief set cluster's level
	* \param [in] level
	*/
	void set_level(size_t level)
	{
		m_level = level;
	}

	/** \brief get cluster's level
	* \return level
	*/
	size_t get_level() const
	{
		return m_level;
	}

	/** \brief set the parent cluster's index
	* \param [in] parent the parent cluster's index
	*/
	void set_parent(size_t parent)
	{
		m_parent = parent;
	}

	/** \brief return parent cluster's index
	 * \return index of parent cluster
	 */
	size_t get_parent() const
	{
		return m_parent;
	}

	/** \brief set the cluster's bounding box
	* \param [in] bb the bounding box
	*/
	void set_bounding_box(bounding_box_t const &bb)
	{
		m_bb = bb;
	}

	/** \brief return cluster's bounding box
	 * \return the bounding box
	 */
	bounding_box_t const &get_bounding_box() const
	{
		return m_bb;
	}

	/** \brief set indices of source nodes
	 * \param [in] node_idx indices of source nodes
	 */
	void set_src_node_idx(idx_list_t const &node_idx)
	{
		m_src_node_idx = node_idx;
	}

	/** \brief set indices of source nodes with range
	 * \param [in] begin begin iterator
	 * \param [in] end end iterator
	 */
	template <class ForwardIt>
	void set_src_node_idx(ForwardIt begin, ForwardIt end)
	{
		m_src_node_idx = idx_list_t(begin, end);
	}

	/** \brief get indices of source nodes
	* \return indices of source nodes
	*/
	idx_list_t const &get_src_node_idx() const
	{
		return m_src_node_idx;
	}

	/** \brief return number of sources contained in cluster
	 * \return the number of source nodes
	 */
	size_t get_n_src_nodes() const
	{
		return m_src_node_idx.size();
	}

	/** \brief set indices of receiver nodes
	* \param [in] node_idx indices of receiver nodes
	*/
	void set_rec_node_idx(idx_list_t const &node_idx)
	{
		m_rec_node_idx = node_idx;
	}

	/** \brief set indices of receiver nodes with range
	 * \param [in] begin begin iterator
	 * \param [in] end end iterator
	 */
	template <class ForwardIt>
	void set_rec_node_idx(ForwardIt begin, ForwardIt end)
	{
		m_rec_node_idx = idx_list_t(begin, end);
	}

	/** \brief get indices of receiver nodes
	 * \return indices of receiver nodes
	 */
	idx_list_t const &get_rec_node_idx() const
	{
		return m_rec_node_idx;
	}

	/** \brief get number of receiver nodes contained in cluster
	* \return number of receiver nodes
	*/
	size_t get_n_rec_nodes() const
	{
		return m_rec_node_idx.size();
	}

	/** \brief add a new child to a cluster
	* \param [in] child the new child to add
	*
	* Does not add child if it is already contained
	*/
	void add_child(size_t child)
	{
		if (std::find(m_children.begin(), m_children.end(), child) ==
			m_children.end())
			m_children.push_back(child);
		else throw std::logic_error("Child already contained");
	}

	/** \brief get indices of children
	* \return indices of children
	*/
	idx_list_t const &get_children() const
	{
		return m_children;
	}

	/** \brief indicates if cluster is leaf
	 * \return true if cluster is leaf
	 * The cluster is leaf if it does not have children
	 */
	bool is_leaf() const
	{
		return m_children.empty();
	}

	/** \brief indicates if cluster is root (0-level)
	 * \return true if cluster is root
	 */
	bool is_root() const
	{
		return m_level == 0;
	}

	/** \brief indicates if cluster is source
	 * \return true if cluster is source
	 * The cluster is source if its source indices is not empty
	 */
	bool is_source() const
	{
		return !m_src_node_idx.empty();
	}

	/** \brief indicates if cluster is receiver
	 * \return true if cluster is receiver
	 * The cluster is receiver if its receiver indices is not empty
	 */
	bool is_receiver() const
	{
		return !m_rec_node_idx.empty();
	}

	/** \brief prints debug information of the cluster to an output stream
	 * \param [in] os the output stream
	 */
	void print_debug(std::ostream &os = std::cout) const
	{
		os << "Level:" << this->get_level() << std::endl;
		os << "Bounding box: " << m_bb << std::endl;
		std::ostream_iterator<size_t> out_it(os, ", ");
		os << "Source indices:" << std::endl;
		std::copy(this->get_src_node_idx().begin(), this->get_src_node_idx().end(), out_it);
		os << std::endl;
		os << "Receiver indices:" << std::endl;
		std::copy(this->get_rec_node_idx().begin(), this->get_rec_node_idx().end(), out_it);
	}

	/** \brief return a cleared (zero) multipole contribution
	 * \return a zero multipole contribution
	 */
	multipole_t zero_multipole() const
	{
		return derived().zero_multipole();
	}

	/** \brief return a cleared (zero) local contribution
	 * \return a zero local contribution
	 */
	local_t zero_local() const
	{
		return derived().zero_local();
	}

private:
	bounding_box_t m_bb;
	size_t m_level;
	size_t m_parent;
	idx_list_t m_children;

	idx_list_t m_src_node_idx;
	idx_list_t m_rec_node_idx;
};

/** \brief inserter operator of a cluster
 * \tparam Derived the CRTP derived class
 * \param [in, out] os the output stream
 * \param [in] c the cluster
 */
template <class Derived>
std::ostream &operator<<(std::ostream &os, cluster_base<Derived> const &c)
{
	c.print_debug(os);
	return os;
}

} // namespace fmm
} // namespace NiHu

#endif // CLUSTER_HPP_INCLUDED
