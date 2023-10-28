#include "spectral_interpolate.h"

namespace NiHu
{
namespace fmm
{

spectral_interpolate::~spectral_interpolate()
{
	destroy();
}

spectral_interpolate::spectral_interpolate(size_t Nto, size_t Nfrom)
	: m_Nto(Nto)
	, m_Nfrom(Nfrom)
	, m_dft_plan(nullptr)
	, m_idft_plan(nullptr)
{
	init();
}

spectral_interpolate::spectral_interpolate(spectral_interpolate const &rhs)
	: m_Nto(0)
	, m_Nfrom(0)
	, m_dft_plan(nullptr)
	, m_idft_plan(nullptr)
{
	*this = rhs;
}

spectral_interpolate::spectral_interpolate(spectral_interpolate &&rhs)
	: m_Nto (rhs.m_Nto)	
	, m_Nfrom (rhs.m_Nfrom)
	, m_dft_plan (rhs.m_dft_plan)
	, m_idft_plan (rhs.m_idft_plan)
	, m_in_spec (std::move(rhs.m_in_spec))
	, m_in (std::move(rhs.m_in))
	, m_out (std::move(rhs.m_out))
	, m_out_spec (std::move(rhs.m_out_spec))
{
	rhs.m_idft_plan = nullptr;
	rhs.m_dft_plan = nullptr;
}

spectral_interpolate &spectral_interpolate::operator=(spectral_interpolate const &rhs)
{
	if (this == &rhs)
		return *this;

	destroy();

	m_Nfrom = rhs.m_Nfrom;
	m_Nto = rhs.m_Nto;

	init();

	return *this;
}

spectral_interpolate &spectral_interpolate::operator=(spectral_interpolate &&rhs)
{
	if (this == &rhs)
		return *this;

	destroy();

	m_Nfrom = rhs.m_Nfrom;
	m_Nto = rhs.m_Nto;

	m_idft_plan = rhs.m_idft_plan;
	m_dft_plan = rhs.m_dft_plan;

	m_in_spec = std::move(rhs.m_in_spec);
	m_in = std::move(rhs.m_in);
	m_out = std::move(rhs.m_out);
	m_out_spec = std::move(rhs.m_out_spec);

	rhs.m_idft_plan = nullptr;
	rhs.m_dft_plan = nullptr;

	return *this;
}

spectral_interpolate::cvector_t const &
spectral_interpolate::interpolate(cvector_t const &from) const
{
	// check argument
	if (from.rows() != 2 * (2 * m_Nfrom) + 1)
		throw std::runtime_error("Invalid size of idft input vector");

	// idft with Nfrom size
	fftw_execute_dft(m_idft_plan, (fftw_complex *)from.data(),
		(fftw_complex *)m_in.data());
	m_in /= std::complex<double>(double(2 * (2 * m_Nfrom) + 1));


	// padding or truncating (m_out has been cleared by init() )
	size_t L = std::min(m_Nfrom, m_Nto);
	m_out.head(L + 1) = m_in.head(L + 1);
	m_out.tail(L) = m_in.tail(L);

	// dft with Nto size
	fftw_execute_dft(m_dft_plan, (fftw_complex *)m_out.data(),
		(fftw_complex *)m_out_spec.data());

	return m_out_spec;
}

void spectral_interpolate::init()
{
	m_in_spec.resize(2 * (2 * m_Nfrom) + 1, 1);
	m_in.resize(2 * (2 * m_Nfrom) + 1, 1);
	m_out.resize(2 * (2 * m_Nto) + 1, 1);
	m_out_spec.resize(2 * (2 * m_Nto) + 1, 1);

	m_idft_plan = fftw_plan_dft_1d(2 * (2 * int(m_Nfrom)) + 1,
		(fftw_complex *)m_in_spec.data(),
		(fftw_complex *)m_in.data(), FFTW_BACKWARD, FFTW_MEASURE);

	m_dft_plan = fftw_plan_dft_1d(2 * (2 * int(m_Nto)) + 1,
		(fftw_complex *)m_out.data(),
		(fftw_complex *)m_out_spec.data(), FFTW_FORWARD, FFTW_MEASURE);

	m_out.setZero();
}

void spectral_interpolate::destroy()
{
	if (m_dft_plan != nullptr)
	{
		fftw_destroy_plan(m_dft_plan);
		m_dft_plan = nullptr;
	}
	if (m_idft_plan != nullptr)
	{
		fftw_destroy_plan(m_idft_plan);
		m_dft_plan = nullptr;
	}
}

} // end of namespace fmm
} // namespace NiHu
