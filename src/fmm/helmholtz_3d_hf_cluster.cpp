/**
 * \file helmholtz_3d_hf_cluster.cpp
 * \brief CLuster implementation for Helmholtz 3D high frequency FMM 
 * \ingroup fmm_helmholtz_3d_hf
 */

#include "helmholtz_3d_hf_cluster.h"

namespace NiHu
{
namespace fmm
{

helmholtz_3d_hf_cluster::multipole_t 
helmholtz_3d_hf_cluster::zero_multipole() const
{
	return multipole_t::Zero(get_level_data().get_unit_sphere().get_s().cols(), 1);
}

helmholtz_3d_hf_cluster::local_t 
helmholtz_3d_hf_cluster::zero_local() const
{
	return local_t::Zero(get_level_data().get_unit_sphere().get_s().cols(), 1);
}

void helmholtz_3d_hf_cluster::set_p_level_data(helmholtz_3d_hf_level_data const * p)
{
	m_p_level_data = p;
}

helmholtz_3d_hf_level_data const &
helmholtz_3d_hf_cluster::get_level_data() const
{
	return *m_p_level_data;
}

} // end of namespace
} // end of namespace NiHu
