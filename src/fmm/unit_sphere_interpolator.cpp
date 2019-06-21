#include "unit_sphere_interpolator.h"

#include <boost/math/special_functions/legendre.hpp>

namespace NiHu
{
namespace fmm
{

interpolator::interpolator(unit_sphere const &Sfrom, unit_sphere const &Sto)
	: m_Lfrom(Sfrom.get_P())
	, m_Lto(Sto.get_P())
	, m_Phi(cmatrix_t::Zero(m_Lfrom, 2 * m_Lfrom))
	, m_Q(cmatrix_t::Zero(m_Lto, 2 * m_Lto))
	, m_res(cvector_t::Zero(m_Lto * 2 * m_Lto, 1))
{
	// plan dft-s
	cmatrix_t a(m_Lfrom, 2 * m_Lfrom), afft(m_Lfrom, 2 * m_Lfrom);
	int n[1];
	n[0] = int(2 * m_Lfrom);
	this->dft_plan = fftw_plan_many_dft(1, n, int(m_Lfrom),
		(fftw_complex *)a.data(), nullptr, int(m_Lfrom), 1,
		(fftw_complex *)afft.data(), nullptr, int(m_Lfrom), 1,
		FFTW_FORWARD, FFTW_MEASURE);

	cmatrix_t bfft(m_Lto, 2 * m_Lto), b(m_Lto, 2 * m_Lto);
	n[0] = int(2 * m_Lto);
	this->idft_plan = fftw_plan_many_dft(1, n, int(m_Lto),
		(fftw_complex *)bfft.data(), nullptr, int(m_Lto), 1,
		(fftw_complex *)b.data(), nullptr, int(m_Lto), 1,
		FFTW_BACKWARD, FFTW_MEASURE);

	// compute A matrices
	int L = int(std::min(m_Lfrom, m_Lto));
	for (int m = -L; m <= L; ++m)
		m_A.push_back(Amatrix(L, m, Sto.get_theta(), Sfrom.get_theta(), Sfrom.get_wtheta()));
}

interpolator::interpolator(interpolator &&other)
	: m_Lfrom(other.m_Lfrom)
	, m_Lto(other.m_Lto)
	, m_A(std::move(other.m_A))
	, dft_plan(other.dft_plan)
	, idft_plan(other.idft_plan)
	, m_Phi(std::move(other.m_Phi))
	, m_Q(std::move(other.m_Q))
	, m_res(std::move(other.m_res))
{
	other.dft_plan = nullptr;
	other.idft_plan = nullptr;
}



interpolator::cvector_t const &
interpolator::interpolate(cvector_t const &other) const
{
	using namespace boost::math::double_constants;

	size_t L = std::min(m_Lfrom, m_Lto);

	// dft along phi
	/** \todo this allocation should be avoided for speed */
	// cmatrix_t Phi(Lfrom, 2 * Lfrom);
	fftw_execute_dft(
		this->dft_plan,
		(fftw_complex *)other.data(),
		(fftw_complex *)m_Phi.data());
	m_Phi *= pi / m_Lfrom;

	// interpolation columnwise with zero padding/truncation
	/** \todo this allocation should be avoided for speed */
	// cmatrix_t Q = cmatrix_t::Zero(Lto, 2 * Lto);
	for (size_t m = 0; m < L; ++m)
		m_Q.col(m) = m_A[m + L] * m_Phi.col(m);
	for (int m = -int(L); m < 0; ++m)
		m_Q.col(2 * m_Lto + m) = m_A[m + L] * m_Phi.col(2 * m_Lfrom + m);

	// idft along phi
	// cvector_t res(Lto * 2 * Lto, 1);
	fftw_execute_dft(
		this->idft_plan,
		(fftw_complex *)m_Q.data(),
		(fftw_complex *)m_res.data());

	return m_res;
}

interpolator const &
interpolator::operator=(interpolator const &other)
{
	if (this == &other)
		return *this;

	if (dft_plan)
		fftw_destroy_plan(dft_plan);
	if (idft_plan)
		fftw_destroy_plan(idft_plan);

	if (other.dft_plan == nullptr)
	{
		dft_plan = nullptr;
		idft_plan = nullptr;
		return *this;
	}

	m_Lfrom = other.m_Lfrom;
	m_Lto = other.m_Lto;
	m_A = other.m_A;
	m_Phi = other.m_Phi;
	m_Q = other.m_Q;
	m_res = other.m_res;

	// plan dft-s
	cmatrix_t a(m_Lfrom, 2 * m_Lfrom), afft(m_Lfrom, 2 * m_Lfrom);
	int n[1];
	n[0] = int(2 * m_Lfrom);
	this->dft_plan = fftw_plan_many_dft(1, n, int(m_Lfrom),
		(fftw_complex *)a.data(), nullptr, int(m_Lfrom), 1,
		(fftw_complex *)afft.data(), nullptr, int(m_Lfrom), 1,
		FFTW_FORWARD, FFTW_MEASURE);
	cmatrix_t bfft(m_Lto, 2 * m_Lto), b(m_Lto, 2 * m_Lto);
	n[0] = int(2 * m_Lto);
	this->idft_plan = fftw_plan_many_dft(1, n, int(m_Lto),
		(fftw_complex *)bfft.data(), nullptr, int(m_Lto), 1,
		(fftw_complex *)b.data(), nullptr, int(m_Lto), 1,
		FFTW_BACKWARD, FFTW_MEASURE);

	return *this;
}

interpolator const &
interpolator::operator=(interpolator &&other)
{
	if (&other == this)
		return *this;

	m_Lfrom = other.m_Lfrom;
	m_Lto = other.m_Lto;
	m_A = std::move(other.m_A);
	dft_plan = other.dft_plan;
	idft_plan = other.idft_plan;
	m_Phi = std::move(other.m_Phi);
	m_Q = std::move(other.m_Q);
	m_res = std::move(other.m_res);

	other.dft_plan = nullptr;
	other.idft_plan = nullptr;

	return *this;
}


interpolator::dmatrix_t interpolator::Amatrix(int L, int m,
	dvector_t const &theta,
	dvector_t const &xi,
	dvector_t const &wxi)
{
	using boost::math::legendre_p;

	dmatrix_t Pi = dmatrix_t::Zero(theta.rows(), L + 1);
	dmatrix_t Pj = dmatrix_t::Zero(xi.rows(), L + 1);
	for (int l = std::abs(m); l <= L; ++l)
	{
		for (int i = 0; i < theta.rows(); ++i)
			Pi(i, l) = legendre_p(l, m, std::cos(theta(i)));
		for (int j = 0; j < xi.rows(); ++j)
			Pj(j, l) = legendre_p(l, m, std::cos(xi(j)));
	}
	dmatrix_t A = dmatrix_t::Zero(theta.rows(), xi.rows());
	for (int l = std::abs(m); l <= L; ++l)
		A += C(l, m) * Pi.col(l) * Pj.col(l).transpose();
	return A * wxi.asDiagonal();
}

double interpolator::C(int l, int m)
{
	double const pi = boost::math::constants::pi<double>();

	double r = (2.0 * l + 1) / (4.0 * pi);
	if (m > 0)
	{
		for (int k = l - m + 1; k <= l + m; ++k)
			r /= k;
	}
	else if (m < 0)
	{
		for (int k = l + m + 1; k <= l - m; ++k)
			r *= k;
	}
	return r;
}

} // end of namespace fmm
} // namespace NiHu