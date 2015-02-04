// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2015  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2015  Peter Rucz <rucz@hit.bme.hu>
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
 * \file elastodynamics_kernel.hpp
 * \ingroup library
 * \brief implementation of fundamental solutions of elastodynamics
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef ELASTODYNAMICS_KERNEL_HPP_INCLUDED
#define ELASTODYNAMICS_KERNEL_HPP_INCLUDED

#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../util/collection.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"
#include "elastostatics_kernel.hpp"

class elastodynamics_data
{
public:
	elastodynamics_data(double nu, double rho, double mu, double omega)
		: m_nu(nu), m_rho(rho), m_mu(mu), m_omega(omega) {}
	double get_poisson_ratio(void) const { return m_nu; }
	double get_mass_density(void) const { return m_rho; }
	double get_shear_modulus(void) const { return m_mu; }
	double get_frequency(void) const { return m_omega; }

private:
	double m_nu;	/**< \brief Poisson's number */
	double m_rho;	/**< \brief mass density */
	double m_mu;	/**< \brief shear modulus */
	double m_omega;	/**< \brief angular frequency */
};

struct DynamicUkernel
{
	typedef Eigen::Matrix<std::complex<double>, 3, 3> return_type;
	
	return_type operator()(
		location_input_3d const &x,
		location_input_3d const &y,
		elastodynamics_data const &data)
	{
		// compute distance and its gradient
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto gradr = rvec.normalized();

		// material properties and Helmholtz numbers
		auto nu = data.get_poisson_ratio();
		auto a = std::sqrt((1.-2.*nu)/2./(1.-nu));
		auto mu = data.get_shear_modulus();
		auto rho = data.get_mass_density();
		auto cS = std::sqrt(mu/rho);

		auto om = data.get_frequency();
		auto kSr = om / cS * r;
		auto kPr = a * kSr;

		// scalar complex helpers
		std::complex<double> const I(0.0, 1.0);

		std::complex<double> psi =
			std::exp(-I*kPr) * (a*a) * (I/kPr + 1./(kPr*kPr)) +
			std::exp(-I*kSr) * (1. - I/kSr - 1./(kSr*kSr));

		std::complex<double> chi =
			std::exp(-I*kPr) * (a*a) * (1. - 3.*I/kPr - 3./(kPr*kPr)) -
			std::exp(-I*kSr) * (1. - 3.*I/kSr - 3./(kSr*kSr));

		// matrix valued result
		return ( psi * return_type::Identity() + chi * (gradr * gradr.transpose()) ) / (4.*M_PI*mu*r);
	}
};

class elastodynamics_3d_U_kernel;

/** \brief the properties of the elastodynamics U kernel */
template <>
struct kernel_traits<elastodynamics_3d_U_kernel>
{
	typedef location_input_3d test_input_t;
	typedef location_input_3d trial_input_t;
	typedef collect<elastodynamics_data> data_t;
	typedef single_brick_wall<DynamicUkernel>::type output_t;
	enum { result_dimension = 3 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	static bool const is_singular = true;
};

/** \brief the singular properties of the elastodynamics U kernel */
template <>
struct singular_kernel_traits<elastodynamics_3d_U_kernel>
{
	typedef asymptotic::inverse<1> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
	typedef elastostatics_3d_U_kernel singular_core_t;
};

class elastodynamics_3d_U_kernel :
	public kernel_base<elastodynamics_3d_U_kernel>
{
public:
	elastodynamics_3d_U_kernel(double nu, double rho, double mu, double omega) :
		kernel_base<elastodynamics_3d_U_kernel>(elastodynamics_data(nu, rho, mu, omega)) {}
};








struct DynamicTkernel
{
	typedef Eigen::Matrix<std::complex<double>, 3, 3> return_type;
	
	return_type operator()(
		location_input_3d const &x,
		location_normal_input_3d const &y,
		elastodynamics_data const &data)
	{
		// compute distance and its gradient
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto const &n = y.get_unit_normal();
		auto gradr = rvec.normalized();
		double rdn = gradr.dot(n);

		// material properties and Helmholtz numbers
		auto nu = data.get_poisson_ratio();
		auto a = std::sqrt((1.-2.*nu)/2./(1.-nu)); /* cs/cp */
		auto a2 = a*a;
		auto mu = data.get_shear_modulus();
		auto lambdapermu = (1./a2 - 2.);
		auto rho = data.get_mass_density();
		auto cS = std::sqrt(mu/rho);

		auto om = data.get_frequency();
		auto kSr = om / cS * r;
		auto kPr = a * kSr;

		// scalar complex helpers
		std::complex<double> const I(0.0, 1.0);

		std::complex<double> psi =
			std::exp(-I*kPr) * a2 * (     I/kPr + 1./(kPr*kPr)) +
			std::exp(-I*kSr) *      (1. - I/kSr - 1./(kSr*kSr));

		std::complex<double> chi =
			std::exp(-I*kPr) * a2 * (1. - 3.*I/kPr - 3./(kPr*kPr)) -
			std::exp(-I*kSr) *      (1. - 3.*I/kSr - 3./(kSr*kSr));

		std::complex<double> dpsi = (
			std::exp(-I*kPr) * a2 * (1.         - 2.*I/kPr - 2./(kPr*kPr)) -
			std::exp(-I*kSr) *      (1. + I*kSr - 2.*I/kSr - 2./(kSr*kSr))
			) / r;

		std::complex<double> dchi = (
			std::exp(-I*kSr) *      (3. + I*kSr - 6.*I/kSr - 6./(kSr*kSr)) -
			std::exp(-I*kPr) * a2 * (3. + I*kPr - 6.*I/kPr - 6./(kPr*kPr))
			) / r;

		std::complex<double> A = dpsi - psi/r;
		std::complex<double> B = dchi - 3.*chi/r;
		std::complex<double> C = chi/r;

		// matrix valued result
		return (
			return_type::Identity() * rdn * (A+C)
			+
			rdn * 2. * B * (gradr * gradr.transpose())
			+
			(2.*C + lambdapermu*(A+B+4.*C)) * (n * gradr.transpose())
			+
			(A+C) * (gradr * n.transpose())
		).transpose() / (4.*M_PI*r);
	}
};

class elastodynamics_3d_T_kernel;

/** \brief the properties of the elastodynamics T kernel */
template <>
struct kernel_traits<elastodynamics_3d_T_kernel>
{
	typedef location_input_3d test_input_t;
	typedef location_normal_input_3d trial_input_t;
	typedef collect<elastodynamics_data> data_t;
	typedef single_brick_wall<DynamicTkernel>::type output_t;
	enum { result_dimension = 3 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = false;
	typedef asymptotic::inverse<2> far_field_behaviour_t;
	static bool const is_singular = true;
};

/** \brief the singular properties of the elastodynamics T kernel */
template <>
struct singular_kernel_traits<elastodynamics_3d_T_kernel>
{
	typedef asymptotic::inverse<2> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
	typedef elastostatics_3d_T_kernel singular_core_t;
};

class elastodynamics_3d_T_kernel :
	public kernel_base<elastodynamics_3d_T_kernel>
{
public:
	elastodynamics_3d_T_kernel(double nu, double rho, double mu, double omega) :
		kernel_base<elastodynamics_3d_T_kernel>(elastodynamics_data(nu, rho, mu, omega)) {}
};



#endif // ELASTODYNAMICS_KERNEL_HPP_INCLUDED

