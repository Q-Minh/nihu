// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file math_functions.hpp
 * \brief general mathematical functions
 */
#ifndef MATH_FUNCTIONS_HPP_INCLUDED
#define MATH_FUNCTIONS_HPP_INCLUDED

#include <cmath>
#include <complex>

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
	/** \brief imaginary unit */
	std::complex<double> const I(0.0, 1.0);
	/** \brief Euler's constant */
	double const gamma(0.57721566490153286060);

	/** \brief small argument expansion of J_nu(z) for nu = 0, 1
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return J_nu(z)
	 */
	template <int nu>
	std::complex<double> J_small(std::complex<double> const &z)
	{
		static_assert(nu == 0 || nu == 1, "unimplemented Bessel J order");

		// upper limit for 1e-8 error
		int N = (int)(3+2.0*std::abs(z));

		std::complex<double> res(1.0), q(z*z/4.0);
		for (int n = N; n > 0; --n)
			res = 1.0 - q*(1.0/(n*(n+nu)))*res;
		if (nu == 1)
			res *= z/2.0;
		return res;
	}

	/** \brief small argument expansion of Y_nu(z)
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return Y_nu(z)
	 */
	template <int nu>
	std::complex<double> Y_small(std::complex<double> const &z);

	template <>
	std::complex<double> Y_small<0>(std::complex<double> const &z)
	{
		int const N = 15;

		std::complex<double> q(z/2.0), q2(q*q);
		std::complex<double> sum(0.0), ss(-1.0);
		double a(0.0);

		for (int k = 1; k <= N; ++k)
		{
			a += 1.0/k;
			ss *= -q2*(1.0/k/k);
			sum += a*ss;
		}

		return 2.0/M_PI*( (std::log(q)+gamma)*J_small<0>(z) + sum );
	}

	template <>
	std::complex<double> Y_small<1>(std::complex<double> const &z)
	{
		int const N = 15;

		std::complex<double> q(z/2.0), q2(q*q);
		std::complex<double> first(2.0*J_small<1>(z)*(std::log(q)+gamma));
		std::complex<double> second(-1.0/q);

		std::complex<double> third(0.0);
		std::complex<double> q2pow(-q);
		double a(1.0);
		double div = 1.0;
		for (int k = 0; k < N; ++k)
		{
			third += div*q2pow*a;

			div /= -(k+1)*(k+2);
			q2pow *= q2;
			a += 1.0/(k+1.0) + 1.0/(k+2.0);
		}

		return 1.0/M_PI * (first + second + third);
	}

	/** \brief large argument expansion of H^(2)_nu(z)
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return H^(2)_nu(z)
	 */
	template <int nu>
	std::complex<double> H_large(std::complex<double> const &z)
	{
		int const N(15);

		std::complex<double> sum(1.0), c(1.0);
		for (int k = 1; k < N; ++k)
		{
			c *= -I/z * ((4.0*nu*nu-(2*k-1)*(2*k-1))/k/8.0);
			sum += c;
		}

		std::complex<double> om(z-nu*M_PI/2.0-M_PI/4.0);
		return std::sqrt(2.0/M_PI/z) * std::exp(-I*om) * sum;
	}

	/** \brief H^(2)_nu(z) Bessel function
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return H^(2)_nu(z)
	 */
	template <int nu>
	std::complex<double> H(std::complex<double> const &z)
	{
		if (std::abs(z) < 6.0)
			return J_small<nu>(z) - I * Y_small<nu>(z);
		else
			return H_large<nu>(z);
	}
}

#endif // MATH_FUNCTIONS_HPP_INCLUDED
