#include "helmholtz_2d_wb_level_data.h"
#include <iostream>
#include <omp.h>

namespace NiHu
{
namespace fmm
{

helmholtz_2d_wb_level_data::helmholtz_2d_wb_level_data()
	: m_interp_ups(1)
	, m_interp_dns(1)
	, m_dft_plan(nullptr)
	, m_multi_paddeds(1)
	, m_multi_spectrums(1)
	, m_idft_plan(nullptr)
	, m_local_paddeds(1)
	, m_locals(1)
	
{
}

helmholtz_2d_wb_level_data::~helmholtz_2d_wb_level_data()
{
	if (m_dft_plan)
		fftw_destroy_plan(m_dft_plan);
	if (m_idft_plan)
		fftw_destroy_plan(m_idft_plan);
}

helmholtz_2d_wb_level_data::helmholtz_2d_wb_level_data(helmholtz_2d_wb_level_data const &rhs)
	: m_high_freq(rhs.m_high_freq)
	, m_interp_ups(rhs.m_interp_ups)
	, m_interp_dns(rhs.m_interp_dns)
	, m_dft_plan(nullptr)
	, m_multi_paddeds(rhs.m_multi_paddeds)
	, m_multi_spectrums(rhs.m_multi_spectrums)
	, m_idft_plan(nullptr)
	, m_local_paddeds(rhs.m_local_paddeds)
	, m_locals(rhs.m_locals)
{
	set_expansion_length(rhs.get_expansion_length());
}


void helmholtz_2d_wb_level_data::init(double drel)
{
	int expansion_length = int(40.0 + 9.0 * drel);
	bool is_high = drel > 3.0;
	set_expansion_length(expansion_length);
	set_high_freq(is_high);
}

void helmholtz_2d_wb_level_data::set_expansion_length(size_t expansion_length)
{
	m_expansion_length = expansion_length;

	// destroy dft's if existing
	if (m_dft_plan)
		fftw_destroy_plan(m_dft_plan);
	if (m_idft_plan)
		fftw_destroy_plan(m_idft_plan);

	size_t S = 2 * (2 * m_expansion_length) + 1;

	// replan dft
	for (auto &a : m_multi_paddeds)
		a.resize(S, 1);
	for (auto &a : m_multi_spectrums)
		a.resize(S, 1);
	m_dft_plan = fftw_plan_dft_1d(int(S), (fftw_complex *)m_multi_paddeds[0].data(),
		(fftw_complex *)m_multi_spectrums[0].data(), FFTW_FORWARD, FFTW_MEASURE);

	// replan idft
	for (auto &a : m_locals)
		a.resize(S, 1);
	for (auto &a : m_local_paddeds)
		a.resize(S, 1);
	m_idft_plan = fftw_plan_dft_1d(int(S), (fftw_complex *)m_locals[0].data(),
		(fftw_complex *)m_local_paddeds[0].data(), FFTW_BACKWARD, FFTW_MEASURE);

	for (auto &a : m_locals)
		a.resize(2 * m_expansion_length + 1, 1);

	for (auto &a : m_multi_paddeds)
		a.setZero();
}

size_t helmholtz_2d_wb_level_data::get_expansion_length() const
{
	return m_expansion_length;
}

size_t helmholtz_2d_wb_level_data::get_data_size() const
{
	size_t M = m_expansion_length;
	if (is_high_freq())
		return 2 * (2 * M) + 1;
	else
		return 2 * M + 1;
}

void helmholtz_2d_wb_level_data::set_high_freq(bool hf)
{
	m_high_freq = hf;
}

bool helmholtz_2d_wb_level_data::is_high_freq() const
{
	return m_high_freq;
}

helmholtz_2d_wb_level_data::cvector_t const &
helmholtz_2d_wb_level_data::dft(cvector_t const &multi) const
{
	auto const idx = omp_get_thread_num();

	if (multi.rows() != 2 * m_expansion_length + 1)
		throw std::runtime_error("Invalid multipole size");

	// pad and shift 
	m_multi_paddeds[idx].tail(m_expansion_length) = multi.head(m_expansion_length);
	m_multi_paddeds[idx].head(m_expansion_length + 1) = multi.tail(m_expansion_length + 1);

	// dft
	fftw_execute_dft(m_dft_plan, (fftw_complex *)m_multi_paddeds[idx].data(),
		(fftw_complex *)m_multi_spectrums[idx].data());

	return m_multi_spectrums[idx];
}

helmholtz_2d_wb_level_data::cvector_t const &
helmholtz_2d_wb_level_data::idft(cvector_t const &local_spect) const
{
	auto const idx = omp_get_thread_num();

	if (local_spect.rows() != 2 * (2 * m_expansion_length) + 1)
		throw std::runtime_error("Invalid local spectrum size");

	// idft
	fftw_execute_dft(m_idft_plan, (fftw_complex *)local_spect.data(),
		(fftw_complex *)m_local_paddeds[idx].data());
	m_local_paddeds[idx] /= double(2 * (2 * m_expansion_length) + 1);

	// shift and depad
	m_locals[idx].head(m_expansion_length) = m_local_paddeds[idx].tail(m_expansion_length);
	m_locals[idx].tail(1 + m_expansion_length) = m_local_paddeds[idx].head(1 + m_expansion_length);

	return m_locals[idx];
}

void helmholtz_2d_wb_level_data::set_interp_up(size_t Nfrom)
{
	m_interp_ups[0] = spectral_interpolate(m_expansion_length, Nfrom);
	for (size_t i = 1; i < m_interp_ups.size(); ++i)
		m_interp_ups[i] = m_interp_ups[0];
}

helmholtz_2d_wb_level_data::cvector_t const &
helmholtz_2d_wb_level_data::interp_up(cvector_t const &multi) const
{
	auto const idx = omp_get_thread_num();
	return m_interp_ups[idx].interpolate(multi);
}

void helmholtz_2d_wb_level_data::set_interp_dn(size_t Nfrom)
{
	m_interp_dns[0] = spectral_interpolate(m_expansion_length, Nfrom);
	for (size_t i = 1; i < m_interp_dns.size(); ++i)
		m_interp_dns[i] = m_interp_dns[0];
}

helmholtz_2d_wb_level_data::cvector_t const &
helmholtz_2d_wb_level_data::interp_dn(cvector_t const &local) const
{
	auto const idx = omp_get_thread_num();
	return m_interp_dns[idx].interpolate(local);
}

void helmholtz_2d_wb_level_data::set_num_threads(size_t n)
{
	// resize all thread_private vectors
	m_interp_dns.resize(n, m_interp_dns[0]);
	m_interp_ups.resize(n, m_interp_ups[0]);
	m_multi_paddeds.resize(n, m_multi_paddeds[0]);
	m_multi_spectrums.resize(n, m_multi_spectrums[0]);
	m_locals.resize(n, m_locals[0]);
	m_local_paddeds.resize(n, m_local_paddeds[0]);
}

} // end of namespace fmm
} // namespace NiHu
