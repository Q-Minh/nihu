/**
 * @file empty_cluster.hpp
 * @brief Definition of the class @ref empty_cluster 
 * @ingroup fmm_clus
 */
#ifndef EMPTY_CLUSTER_HPP_INCLUDED
#define EMPTY_CLUSTER_HPP_INCLUDED

#include "cluster.hpp"

namespace NiHu
{
namespace fmm
{

template <size_t Dim>
class empty_cluster;

/** 
 * @brief Traits of empty cluster 
 * @tparam Dim Space dimension
 */
template <size_t Dim>
struct cluster_traits<empty_cluster<Dim> >
{
	static size_t const dimension = Dim;
	typedef void multipole_t;
	typedef void local_t;
};

/**
 * @brief Empty cluster class 
 * @tparam[in] Dim Space dimension
 */
template <size_t Dim>
class empty_cluster
	: public cluster_base<empty_cluster<Dim> >
{
public:
	typedef cluster_base<empty_cluster<Dim> > base_t;
	typedef typename base_t::multipole_t multipole_t;
	typedef typename base_t::local_t local_t;

	multipole_t zero_multipole() const
	{
		return;
	}

	local_t zero_local() const
	{
		return;
	}
};

} // end of namespace fmm
} // end of namespace NiHu


#endif /* EMPTY_CLUSTER_HPP_INCLUDED */
