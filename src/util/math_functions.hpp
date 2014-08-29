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

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <complex>
#include <stdexcept>
#include <iostream>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

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
		return 1. - x*x/6. * (1. - x*x/20.);
}

/** \brief namespace encapsulating Bessel functions */
namespace bessel
{
	/** \brief metafunction converting a floating point type to a complex type */
	template <class T>
	struct make_complex { typedef std::complex<T> type; };

	template <class T>
	struct make_complex<std::complex<T> > { typedef std::complex<T> type; };

	/** \brief limit between small and large argument series */
	double const large_lim(7.);
	/** \brief imaginary unit */
	std::complex<double> const I(0., 1.);
	/** \brief Euler's constant */
	double const gamma(0.57721566490153286060);

	/** \brief large argument Bessel function Taylor series coefficients
	 * \tparam T the argument type
	 * \param [in] nu the Bessel function order
	 * \param [in] z the Bessel argument
	 * \param [out] u the magnitude approximation
	 * \param [out] u the phase approximation
	 */
	template <class T>
	void mag_arg_large(int nu, T const &z, T &u, T &phi)
	{
		double mag[2][7] = {
			{
				1.,
				-6.25e-2,
				1.03515625e-1,
				-5.428466796875e-1,
				5.8486995697021484375,
				-1.06886793971061706543e2,
				2.96814293784275650978e3
			},
			{
				1.,
				1.875e-1,
				-1.93359375e-1,
				8.052978515625e-1,
				-7.7399539947509765625,
				1.32761824250221252441e2,
				-3.54330366536602377892e3
			}
		};
		double arg[2][7] = {
			{
				-1.25e-1,
				6.51041666666666712926e-2,
				-2.09570312499999994449e-1,	
				1.63806588309151779370,
				-2.34751277499728736586e1,
				5.35640519510615945364e2,
				-1.78372796889474739146e4
			},
			{
				3.75e-1,
				-1.640625e-1,
				3.70898437500000011102e-1,
				-2.36939784458705338110,
				3.06240119934082031250e1,
				-6.59185221823778988437e2,
				2.11563140455278044101e4
			}
		};

		u = phi = T(0.);
		auto q(1./z/z);
		for (int m = sizeof(mag[0])/sizeof(double)-1; m >= 0; --m)
		{
			u = q * u + mag[nu][m];
			phi = q * phi + arg[nu][m];
		}
		phi = phi/z + z - (nu/2.+.25)*M_PI;
	}

	/** \brief small argument expansion of J_nu(z) for nu = 0, 1
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return J_nu(z)
	 */
	template <int nu, class T>
	T J_small(T const &z)
	{
		static_assert(nu == 0 || nu == 1, "unimplemented Bessel J order");

		// upper limit for 1e-8 error
		int N = (int)(3+2.*std::abs(z));

		T res(1.), q(z*z/4.);
		for (int n = N; n > 0; --n)
			res = 1. - q*(1./(n*(n+nu)))*res;
		if (nu == 1)
			res *= z/2.;
		return res;
	}

	/** \brief large argument expansion of J_nu(z) for nu = 0, 1
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return J_nu(z)
	 */
	template <int nu, class T>
	T J_large(T const &z)
	{
		T mag, arg;
		if (std::real(z) < 0)
		{
			double const C(nu==0 ? 1. : -1);
			mag_arg_large(nu, -z, mag, arg);
			return std::sqrt(2./(M_PI * -z)) * mag * std::cos(arg) * C;
		}
		else
		{
			mag_arg_large(nu, z, mag, arg);
			return std::sqrt(2./(M_PI * z)) * mag * std::cos(arg);
		}
	}

	/** \brief Bessel function J_nu(z) for nu = 0, 1
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return J_nu(z)
	 */
	template <int nu, class T>
	T J(T const &z)
	{
		if (std::abs(z) < large_lim)
			return J_small<nu>(z);
		else
			return J_large<nu>(z);
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

	/** \brief large argument expansion of Y_nu(z) for nu = 0, 1
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return Y_nu(z)
	 */
	template <int nu>
	std::complex<double> Y_large(std::complex<double> const &z)
	{
		std::complex<double> mag, arg;
		if (std::real(z) < 0)
		{
			double const C1(nu == 0 ? 1. : -1.);
			std::complex<double> const C2(0., std::imag(z) < 0 ? -2. : 2.);
			mag_arg_large(nu, -z, mag, arg);
			return std::sqrt(2./(M_PI * -z)) * mag * (std::sin(arg) + C2 * std::cos(arg)) * C1;
		}
		else
		{
			mag_arg_large(nu, z, mag, arg);
			return std::sqrt(2./(M_PI * z)) * mag * std::sin(arg);
		}
	}


	/** \brief Bessel function Y_nu(z) for nu = 0, 1
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return Y_nu(z)
	 */
	template <int nu>
	std::complex<double> Y(std::complex<double> const &z)
	{
		if (std::abs(z) < large_lim)
			return Y_small<nu>(z);
		else
			return Y_large<nu>(z);
	}

	/** \brief large argument expansion of H^(2)_nu(z)
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return H^(2)_nu(z)
	 */
	template <int nu, int kind, class T>
	typename make_complex<T>::type H_large(T const &z)
	{
		static_assert((kind == 1) || (kind == 2), "invalid kind argument of bessel::H");
		
		if (std::real(z) < 0 && std::abs(std::imag(z)) < 8.)
			throw std::runtime_error("NiHu");

		double const C = (kind == 2 ? -1. : 1.);

		typename make_complex<T>::type mag, arg;
		mag_arg_large(nu, z, mag, arg);
		
		return std::sqrt(2./M_PI/z) * mag * std::exp(C*I*arg);
	}

	/** \brief H^(2)_nu(z) Bessel function
	 * \tparam nu the Bessel function's order
	 * \tparam kind the Bessel function's kind
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return H^(2)_nu(z)
	 */
	template <int nu, int kind, class T>
	typename make_complex<T>::type H(T const &z)
	{
		static_assert(kind == 1 || kind == 2, "invalid kind argument of bessel::H");
		double const C = (kind == 2 ? -1. : 1.);
		if (std::abs(z) < large_lim)
			return J_small<nu>(z) + C * I * Y_small<nu>(z);
		else
			return H_large<nu, kind>(z);
	}

	/** \brief K_nu(z) modified Bessel function
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return K_nu(z)
	 */
	template <int nu, class T>
	typename make_complex<T>::type K(T const &z);

	template <>
	std::complex<double> K<0>(std::complex<double> const &z)
	{
		if (std::imag(z) < 0.)
			return .5*I*M_PI*H<0, 1>(I*z);
		else
			return -.5*I*M_PI*H<0, 2>(-I*z);
	}

	template <>
	std::complex<double> K<1>(std::complex<double> const &z)
	{
		if (std::imag(z) < 0.)
			return -.5*M_PI*H<1, 1>(I*z);
		else
			return -.5*M_PI*H<1, 2>(-I*z);
	}
}

#endif // MATH_FUNCTIONS_HPP_INCLUDED

