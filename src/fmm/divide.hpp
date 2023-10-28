/** 
 * @file divide.hpp
 * @brief Cluster division strategies
 * @ingroup fmm_clus
 */

#ifndef DIVIDE_HPP_INCLUDED
#define DIVIDE_HPP_INCLUDED

#include "cluster.hpp"
#include "util/crtp_base.hpp"

#include <algorithm> // std::max

namespace NiHu
{
namespace fmm
{

/**
 * @brief Base CRTP class for cluster division
 * @details
 * Cluster division classes are predicates derived from this base class using 
 * static polymorphism. Each derived class should implement the template
 * operator() method.
 */
template <class Derived>
class divide_base
{
public:
	NIHU_CRTP_HELPERS
	
	template <class Cluster>
	bool operator()(Cluster const &c) const
	{
		return derived().operator()(c);
	}
};
	
/** @brief class representing a balanced tree division predicate */
class divide_depth : public divide_base<divide_depth>
{
	/** 
	 * @brief Maximal division depth
	 * @details 
	 * The resulting number of levels is depth + 1
	 */
	size_t m_depth;

public:
	/** 
	 * @brief Constructor
	 * @param [in] depth Maximal depth of the tree
	 */
	divide_depth(size_t depth)
		: m_depth(depth)
	{
	}

	/** 
	 * @brief Determine if a cluster needs to be divided
	 * @tparam Cluster Cluster type
	 * @param [in] c the cluster to divide
	 * @return @c true if the cluster needs to be divided
	 */
	template <class Cluster>
	bool operator()(Cluster const &c) const
	{
		return c.get_level() < m_depth;
	}
};

/** \brief class representing a cluster division based on number of nodes */
class divide_num_nodes : public divide_base<divide_num_nodes>
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
	bool operator()(Cluster const &c) const
	{
		size_t s = c.get_n_src_nodes();
		size_t r = c.get_n_rec_nodes();
		return std::max(s, r) > m_max_nodes;
	}
};

/** \brief class representing a balanced tree division predicate by leaf diameter */
class divide_diameter : public divide_base<divide_diameter>
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
	bool operator()(Cluster const &c) const
	{
		return c.get_bounding_box().get_diameter() > m_diameter;
	}
};

} // end of namespace fmm
} // end of namespace NiHu

#endif /* DIVIDE_H_INCLUDED */
