/// \file helmholtz_3d_hf_cluster.h
/// \brief declaration of clas fmm::helmholtz_3d_hf_cluster

#ifndef FMM_HELMHOLTZ_3D_HF_CLUSTER_H_INCLUDED
#define FMM_HELMHOLTZ_3D_HF_CLUSTER_H_INCLUDED

#include "cluster.hpp"
#include "helmholtz_3d_hf_level_data.h"

namespace NiHu
{
namespace fmm
{

/// \brief cluster type of the Helmholtz 3D High frequency FMM
class helmholtz_3d_hf_cluster;

/// \brief specialisation of cluster traits for the Helmholtz 3D HF cluster
template <>
struct cluster_traits<helmholtz_3d_hf_cluster>
{
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
	static size_t const dimension = 3;
	typedef cvector_t multipole_t;
	typedef cvector_t local_t;
};

/// \brief cluster type of the Helmholtz 3D High frequency FMM
class helmholtz_3d_hf_cluster
	: public cluster_base<helmholtz_3d_hf_cluster>
{
	typedef cluster_base<helmholtz_3d_hf_cluster> base_t;

public:
	/// \brief the space dimension
	static size_t const dimension = base_t::dimension;
	/// \brief the multipole type
	typedef typename base_t::multipole_t multipole_t;
	/// \brief the local type
	typedef typename base_t::local_t local_t;

	/// \brief return a zero (cleared) multipole contribution
	multipole_t zero_multipole() const;

	/// \brief return a zero (cleared) local contribution
	local_t zero_local() const;

	/// \brief set the pointer to the cluster's level data
	/// \param [in] p pointer to the level data
	void set_p_level_data(helmholtz_3d_hf_level_data const * p);

	/// \brief return the cluster's level data
	/// \return the cluster's level data
	helmholtz_3d_hf_level_data const &get_level_data() const;

private:
	helmholtz_3d_hf_level_data const *m_p_level_data;
};

} // end of namespace fmm
} // namespace NiHu

#endif //FMM_HELMHOLTZ_3D_HF_CLUSTER_H_INCLUDED
