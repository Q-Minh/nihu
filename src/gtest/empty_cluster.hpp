#ifndef EMPTY_CLUSTER_HPP_INCLUDED
#define EMPTY_CLUSTER_HPP_INCLUDED

#include "fmm/cluster.hpp"

template <size_t Dim>
class empty_cluster;

namespace NiHu
{
namespace fmm
{
template <size_t Dim>
struct cluster_traits<empty_cluster<Dim> >
{
	static size_t const dimension = Dim;
	typedef void multipole_t;
	typedef void local_t;
};
}
}


template <size_t Dim>
class empty_cluster
	: public NiHu::fmm::cluster_base<empty_cluster<Dim> >
{
public:
	typedef NiHu::fmm::cluster_base<empty_cluster<Dim> > base_t;
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

#endif
