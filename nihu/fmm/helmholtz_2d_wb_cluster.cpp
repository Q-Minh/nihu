/**
 * \file helmholtz_2d_wb_cluster.cpp
 * \brief Cluster for the 2D Wideband Helmholtz FMM
 * \ingroup fmm_helmholtz_2d_wb
 */

#include "helmholtz_2d_wb_cluster.h"

namespace NiHu
{
namespace fmm
{

void helmholtz_2d_wb_cluster::set_p_level_data(helmholtz_2d_wb_level_data const *p)
{
	m_p_level_data = p;
}

helmholtz_2d_wb_level_data const &
helmholtz_2d_wb_cluster::get_level_data() const
{
	return *m_p_level_data;
}

helmholtz_2d_wb_cluster::multipole_t 
helmholtz_2d_wb_cluster::zero_multipole() const
{
	return multipole_t::Zero(get_level_data().get_data_size(), 1);
}

helmholtz_2d_wb_cluster::local_t 
helmholtz_2d_wb_cluster::zero_local() const
{
	return local_t::Zero(get_level_data().get_data_size(), 1);
}

} // end of namespace fmm
} // end of namespace NiHu
