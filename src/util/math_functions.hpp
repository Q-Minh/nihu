/**
 * \file math_functions.hpp
 * \brief general mathematical functions
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef MATH_FUNCTIONS_HPP_INCLUDED
#define MATH_FUNCTIONS_HPP_INCLUDED

/**
 * \brief \f$ sinc(x) = \sin(x) / x \f$ function
 * \tparam T type of x
 * \param [in] x the input
 * \return sinc(x)
 */
template <class T>
static T sinc(T const &x)
{
	if (std::abs(x) > 1e-3)
		return std::sin(x) / x;
	else
		return 1.0 - x*x/6.0 * (1.0 - x*x/20.0);
}

namespace bessel
{
	std::complex<double> const I(0.0, 1.0);

	template <int nu>
	std::complex<double> om_func(std::complex<double> const &z)
	{
		return z-nu*M_PI/2.0-M_PI/4.0;
	}

	template <int nu>
	std::complex<double> J_small(std::complex<double> const &z);

	template <>
	std::complex<double> J_small<0>(std::complex<double> const &z)
	{
		std::complex<double> res(1.0), q(z*z/4.0);
		for (int n = 10; n > 0; --n)
			res = 1.0 - q*(1.0/(n*n))*res;
		return res;
	}

	template <int nu>
	std::complex<double> Y_small(std::complex<double> const &z);

	template <>
	std::complex<double> Y_small<0>(std::complex<double> const &z)
	{
		int const N = 10;
		double const gamma(0.57721566490153286060);

		std::complex<double> q(-z*z/4.0);
		std::complex<double> sum(0.0), ss(-1.0);
		double a(0.0);

		for (int k = 1; k <= N; ++k)
		{
			a += 1.0/k;
			ss *= q*(1.0/k/k);
			sum += a*ss;
		}

		return 2.0/M_PI*(
			(std::log(z/2.0)+gamma)*J_small<0>(z) + sum
		);
	}

	template <int nu>
	std::complex<double> H_large(std::complex<double> const &z)
	{
		int const N(2);
		std::complex<double> s(0.0), c(1.0);
		double a(1.0);
		for (int k = 0; k < N; ++k)
		{
			if (k > 0)
				a *= (4*nu*nu-(2*k-1)*(2*k-1))/k/8.0;
			s += c * a;
			c *= I/z;
		}

		return std::sqrt(2.0/M_PI/z) *
			std::exp(I*om_func<nu>(z)) *
			s;
	}

	template <int nu>
	std::complex<double> H(std::complex<double> const &z)
	{
		if (std::abs(z) < 5.0)
			return J_small<nu>(z) + I * Y_small<nu>(z);
		else
			return H_large<nu>(z);
	}
}

#endif // MATH_FUNCTIONS_HPP_INCLUDED
