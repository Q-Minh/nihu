/** \file divide.h
 * \brief cluster division strategies
 */
#ifndef DIVIDE_H_INCLUDED
#define DIVIDE_H_INCLUDED

#include "cluster.hpp"

#include <algorithm> // std::max

namespace NiHu
{
namespace fmm
{

/** \brief class representing a balanced tree division predicate */
class divide_depth
{
	size_t m_depth;

public:
	/** \brief constructor
	 * \param [in] depth maximal depth of the tree
	 */
	divide_depth(size_t depth)
		: m_depth(depth)
	{
	}

	/** \brief determine if a cluster needs to be divided or not
	 * \tparam Cluster the cluster type
	 * \param [in] c the cluster to divide
	 * \return true if the cluster needs to be divided
	 */
	template <class Cluster>
	bool operator()(Cluster const &c)
	{
		return c.get_level() < m_depth;
	}
};

/** \brief class representing a cluster division based on number of nodes */
class divide_num_nodes
{
	size_t m_max_nodes;

public:
	/** \brief constructor
	 * \param [in] max_nodes maximal number of nodes in a leaf cluster
	 */
	divide_num_nodes(size_t max_nodes)
		: m_max_nodes(max_nodes)
	{
	}

	/** \brief determine if a cluster needs to be divided
	 * \tparam Cluster the cluster type
	 * \param [in] c the cluster to divide
	 * \return true if the cluster needs to be divided
	 */
	template <class Cluster>
	bool operator()(Cluster const &c)
	{
		size_t s = c.get_n_src_nodes();
		size_t r = c.get_n_rec_nodes();
		return std::max(s, r) > m_max_nodes;
	}
};

/** \brief class representing a balanced tree division predicate by leaf diameter */
class divide_diameter
{
	double m_diameter;

public:
	/** \brief constructor
	 * \param [in] d diameter limit
	 */
	divide_diameter(double d)
		: m_diameter(d)
	{
	}

	/** \brief determine if a cluster needs to be divided or not
	 * \tparam Cluster the cluster type
	 * \param [in] c the cluster to divide
	 * \return true if the cluster needs to be divided
	 */
	template <class Cluster>
	bool operator()(Cluster const &c)
	{
		return c.get_bounding_box().get_diameter() > m_diameter;
	}
};

} // end of namespace fmm
} // namespace NiHu

#endif // DIVIDE_H_INCLUDED
