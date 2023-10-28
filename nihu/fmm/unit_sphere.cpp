#include "unit_sphere.h"

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/constants/constants.hpp>

#include "nihu/core/gaussian_quadrature.hpp"


namespace NiHu
{
namespace fmm
{

unit_sphere::unit_sphere(size_t P)
	: m_P(P)
{
	using namespace boost::math::double_constants;
	auto res = NiHu::gauss_impl<double>(unsigned(m_P));
	auto const &xi = res.col(0);
	this->wtheta = res.col(1);

	// number of quadrature points
	size_t N = 2 * m_P*m_P;
	s.resize(3, N);
	w.resize(N, 1);

	theta = acos(xi.array());
	phi.resize(2 * m_P);

	for (size_t i = 0; i < 2 * m_P; ++i)
	{
		phi(i) = i * pi / m_P;
		s.block(0, i*m_P, 1, m_P) = std::cos(phi(i)) * sin(theta.array().transpose());
		s.block(1, i*m_P, 1, m_P) = std::sin(phi(i)) * sin(theta.array().transpose());
		s.block(2, i*m_P, 1, m_P) = xi.transpose();
		w.segment(i*m_P, m_P) = this->wtheta * pi / m_P;
	}
}

unit_sphere::cmatrix_t unit_sphere::spht(cvector_t const &f, int L) const
{
	using boost::math::spherical_harmonic;
	Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Flm(L + 1, 2 * L + 1);
	Flm.setZero();
	for (int l = 0; l <= L; ++l)
		for (int m = -l; m <= l; ++m)
			for (int k = 0; k < s.cols(); ++k)	// integration loop
				Flm(l, m + L) +=
				conj(spherical_harmonic(l, m, theta(k%m_P), phi(k / m_P)))
				* f(k) * w(k);
	return Flm;
}

unit_sphere::cvector_t unit_sphere::ispht(cmatrix_t const &Flm) const
{
	using boost::math::spherical_harmonic;
	Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> f(s.cols(), 1);
	f.setZero();
	int L = int(Flm.rows() - 1);
	for (int l = 0; l <= L; ++l)
		for (int m = -l; m <= l; ++m)
			for (int k = 0; k < s.cols(); ++k) // integration loop
				f(k) += Flm(l, m + L)
				* spherical_harmonic(l, m, theta(k%m_P), phi(k / m_P));
	return f;
}

} // namespace fmm
} // namespace NiHu