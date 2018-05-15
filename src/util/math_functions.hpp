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

#include "math_constants.hpp"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <iostream>

namespace NiHu
{

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

	/** \brief large argument Bessel function Taylor series coefficients
	 * \tparam T the argument type
	 * \param [in] nu the Bessel function order
	 * \param [in] z the Bessel argument
	 * \param [out] u the magnitude approximation
	 * \param [out] phi the phase approximation
	 */
	template <class T>
	void mag_arg_large(int nu, T const &z, T &u, T &phi)
	{
		double mag[3][7] = {
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
			},
			{
				1,
				9.3750000000e-01,
				7.9101562500e-01,
				-3.0487060547e+00,
				1.9199895859e+01,
				-2.5916166008e+02,
				6.0841126316e+03
			}
		};
		double arg[3][7] = {
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
			},
			{
				1.8750000000e0,
				-3.5156250000e-1,
				-1.4501953125e0,
				8.3074297224e0,
				-7.1739864349e1,
				1.2458673851e3,
				-3.5571449717e4				
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

	/** \brief small argument expansion of J_nu(z) for nu = 0, 1, 2
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return J_nu(z)
	 */
	template <int nu, class T>
	T J_small(T const &z)
	{
		// upper limit for 1e-8 error
		int N = (int)(3+2.*std::abs(z));

		T q = z/2.;

		T res(1.), q2(q*q);
		for (int k = N; k > 0; --k)
			res = 1. - q2*(1./(k*(k+nu)))*res;
		for (int k = 1; k <= nu; ++k)
			res *= q/T(k);
		return res;
	}

	/** \brief large argument expansion of J_nu(z) for nu = 0, 1, 2
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return J_nu(z)
	 */
	template <int nu, class T>
	T J_large(T const &z)
	{
		static_assert(nu >= 0 && nu <= 2 , "unimplemented Bessel J order");

		T mag, arg;
		if (std::real(z) < 0)
		{
			double const C(nu%2==0 ? 1. : -1);
			mag_arg_large(nu, -z, mag, arg);
			return std::sqrt(2./(M_PI * -z)) * mag * std::cos(arg) * C;
		}
		else
		{
			mag_arg_large(nu, z, mag, arg);
			return std::sqrt(2./(M_PI * z)) * mag * std::cos(arg);
		}
	}

	/** \brief Bessel function J_nu(z)
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return J_nu(z)
	 */
	template <int nu, class T>
	T J(T const &z)
	{
		 return std::abs(z) < large_lim ? J_small<nu>(z) : J_large<nu>(z);
	}


	/** \brief small argument expansion of Y_nu(z)
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return Y_nu(z)
	 */
	template <int nu>
	std::complex<double> Y_small(std::complex<double> const &z)
	{
		// upper limit for 1e-8 error
		int const N = (int)(4+2.*std::abs(z));

		std::complex<double> q(z/2.0), q2(q*q);
		std::complex<double> first(2.0*J_small<nu>(z)*(std::log(q)+M_EULER_GAMMA));
		std::complex<double> second;
		switch (nu)
		{
			case 0: second = 0.; break;
			case 1: second = -1./q; break;
			case 2: second = -1. -1./q2; break;
		}

		std::complex<double> third(0.0);
		std::complex<double> q2pow(-std::pow(q, nu));
		
		double a(0.);
		for (int k = 1; k <= nu; ++k)
			a += 1./k;
		
		double div = 1.;
		for (int k = 1; k <= nu; ++k)
			div /= k;
		
		for (int k = 0; k < N; ++k)
		{
			third += div*q2pow*a;

			div /= -(k+1)*(k+nu+1);
			q2pow *= q2;
			a += 1./(k+1.) + 1./(k+nu+1.);
		}

		return 1./M_PI * (first + second + third);
	}

	/** \brief large argument expansion of Y_nu(z) for nu = 0, 1, 2
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return Y_nu(z)
	 */
	template <int nu>
	std::complex<double> Y_large(std::complex<double> const &z)
	{
		static_assert(nu >= 0 || nu <= 2, "unimplemented Bessel Y order");
		
		std::complex<double> mag, arg;
		if (std::real(z) < 0)
		{
			double const C1(nu%2 == 0 ? 1. : -1.);
			std::complex<double> const C2(0., std::imag(z) < 0 ? -2. : 2.);
			std::complex<double> const MAG(std::sqrt(C2*C2+1.));
			std::complex<double> const ARG(std::atan(1./C2));
			mag_arg_large(nu, -z, mag, arg);
			return std::sqrt(2./(M_PI * -z)) * mag * MAG * std::cos(arg-ARG) * C1;
		}
		else
		{
			mag_arg_large(nu, z, mag, arg);
			return std::sqrt(2./(M_PI * z)) * mag * std::sin(arg);
		}
	}


	/** \brief Bessel function Y_nu(z)
	 * \tparam nu the Bessel function's order
	 * \param [in] z the argument
	 * \return Y_nu(z)
	 */
	template <int nu>
	std::complex<double> Y(std::complex<double> const &z)
	{
		return std::abs(z) < large_lim ? Y_small<nu>(z) : Y_large<nu>(z);
	}

	/** \brief large argument expansion of H^(K)_nu(z)
	 * \tparam nu the Bessel function's order
	 * \tparam T the Bessel argument type
	 * \param [in] z the argument
	 * \return H^(2)_nu(z)
	 */
	template <int nu, int kind, class T>
	typename make_complex<T>::type H_large(T const &z)
	{
		static_assert((kind == 1) || (kind == 2), "invalid kind argument of bessel::H");
		double const C = (kind == 2 ? -1. : 1.);

		typename make_complex<T>::type mag, arg;

		if (std::real(z) < 0 && std::abs(std::imag(z)) < 8.)
		{
			double const C1(nu%2 == 0 ? 1. : -1.);

			mag_arg_large(nu, -z, mag, arg);

			double sgn = std::imag(z) < 0 ? -1. : 1.;
			double MAG(-sgn*std::sqrt(3.));
			std::complex<double> const ARG(std::atan(-sgn*I*.5));

			return std::sqrt(2./(M_PI * -z)) * mag * (std::cos(arg) + C * MAG * std::cos(arg-ARG)) * C1;
		}

		mag_arg_large(nu, z, mag, arg);
		
		return std::sqrt(2./M_PI/z) * mag * std::exp(C*I*arg);
	}

	/** \brief H^(K)_nu(z) Bessel function
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
}

}

#endif // MATH_FUNCTIONS_HPP_INCLUDED

