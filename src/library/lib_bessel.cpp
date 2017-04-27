#include "../util/math_functions.hpp"


namespace NiHu
{
namespace bessel
{
	template <>
	std::complex<double> Y_small<0>(std::complex<double> const &z)
	{
		// upper limit for 1e-8 error
		int const N = (int)(4+2.*std::abs(z));

		std::complex<double> q(z/2.), q2(q*q);
		std::complex<double> sum(0.), ss(-1.);
		double a(0.);

		for (int k = 1; k <= N; ++k)
		{
			a += 1./k;
			ss *= -q2*(1./k/k);
			sum += a*ss;
		}

		return 2.0/M_PI*( (std::log(q)+gamma)*J_small<0>(z) + sum );
	}

	template <>
	std::complex<double> Y_small<1>(std::complex<double> const &z)
	{
		// upper limit for 1e-8 error
		int const N = (int)(4+2.*std::abs(z));

		std::complex<double> q(z/2.0), q2(q*q);
		std::complex<double> first(2.0*J_small<1>(z)*(std::log(q)+gamma));
		std::complex<double> second(-1.0/q);

		std::complex<double> third(0.0);
		std::complex<double> q2pow(-q);
		double a(1.);
		double div = 1.;
		for (int k = 0; k < N; ++k)
		{
			third += div*q2pow*a;

			div /= -(k+1)*(k+2);
			q2pow *= q2;
			a += 1./(k+1.) + 1./(k+2.);
		}

		return 1./M_PI * (first + second + third);
	}

	template <>
	std::complex<double> K<0>(std::complex<double> const &z)
	{
		return .5*I*M_PI* (std::imag(z) < 0. ? H<0, 1>(I*z) : -H<0, 2>(-I*z));
	}

	template <>
	std::complex<double> K<1>(std::complex<double> const &z)
	{
		return -.5*M_PI * (std::imag(z) < 0. ? H<1, 1>(I*z) : H<1, 2>(-I*z));
	}
}
}
