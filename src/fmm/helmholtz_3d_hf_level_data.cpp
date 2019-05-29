#include "helmholtz_3d_hf_level_data.h"

#include <omp.h>


namespace NiHu
{
namespace fmm
{

void helmholtz_3d_hf_level_data::set_expansion_length(size_t L)
{
	m_expansion_length = L;
	m_unit_sphere = unit_sphere(L);
}

unit_sphere const &
helmholtz_3d_hf_level_data::get_unit_sphere() const
{
	return m_unit_sphere;
}

size_t helmholtz_3d_hf_level_data::get_expansion_length() const
{
	return m_expansion_length;
}

void helmholtz_3d_hf_level_data::set_num_threads(size_t n)
{
	m_interp_dns.resize(n, m_interp_dns[0]);
	m_interp_ups.resize(n, m_interp_ups[0]);
}

helmholtz_3d_hf_level_data::cvector_t const &
helmholtz_3d_hf_level_data::interp_up(cvector_t const &x) const
{
	auto const idx = omp_get_thread_num();
	return m_interp_ups[idx].interpolate(x);
}

helmholtz_3d_hf_level_data::cvector_t const &
helmholtz_3d_hf_level_data::interp_down(cvector_t const &x) const
{
	auto const idx = omp_get_thread_num();
	return m_interp_dns[idx].interpolate(x);
}

void helmholtz_3d_hf_level_data::set_interp_up(interpolator const &interp)
{
	for (auto &a : m_interp_ups)
		a = interp;
}

void helmholtz_3d_hf_level_data::set_interp_dn(interpolator const &interp)
{
	for (auto &a : m_interp_dns)
		a = interp;
}


} // end of namespace fmm
} // namespace NiHu