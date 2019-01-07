// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2018  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2018  Peter Rucz <rucz@hit.bme.hu>
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

/** \file helmholtz_singular_integrals.hpp
 * \brief (Semi)analytical expressions for the singular integrals of Helmholtz kernels
 */
#ifndef NIHU_HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
#define NIHU_HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "helmholtz_kernel.hpp"
#include "guiggiani_1992.hpp"
#include "lib_element.hpp"
#include "laplace_singular_integrals.hpp"
#include "normal_derivative_singular_integrals.hpp"
#include "plane_element_helper.hpp"
#include "../core/match_types.hpp"
#include "../core/singular_integral_shortcut.hpp"
#include "../util/math_functions.hpp"

#if NIHU_MEX_DEBUGGING
#include <mex.h>
#endif

namespace NiHu
{

/** \brief Collocational integral of the 2D SLP kernel over a general curved line with general shape sets
 * Singularity subtraction technique is used.
 * The singularty is subtracted in the reference domain, and is integrated analytically.
 * The remaining regular part is integrated numerically with standard Gaussian quadratures
 */
template <class TestField, class TrialField, size_t order>
class helmholtz_2d_SLP_collocation_general
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	
	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;
	
	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;
	
	typedef Eigen::Matrix<std::complex<double>, nTest, nTrial> result_t;
	
	typedef typename trial_field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;
	
	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	
	typedef regular_quad_store<domain_t, order> quadr_t;
	
public:
	template <class WaveNumber>
	static result_t eval(elem_t const &elem, WaveNumber const &k)
	{
		result_t result = result_t::Zero();
		
		helmholtz_2d_SLP_kernel<WaveNumber> kernel(k);
		
		xi_t const &a = domain_t::get_corner(0);
		xi_t const &b = domain_t::get_corner(1);
		
		// traverse collocation points
		for (size_t i = 0; i < nTest; ++i)
		{
			xi_t const &xi0 = test_shape_t::corner_at(i);
			
			// trial shape function at the singular point
			auto N0 = trial_shape_t::template eval_shape<0>(xi0);
			
			// singular point and normal at singular point
			x_t x = elem.get_x(xi0);
			x_t Jxvec = elem.get_normal(xi0);
			double jac0 = Jxvec.norm();
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [a xi0]
				xi_t xi = it->get_xi() * ((xi0(0)-a(0))/2.) + (xi0+a)/2.;
				double w = it->get_w() * ((xi0(0)-a(0))/2.);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y);
				// Shape function
				auto N = trial_shape_t::template eval_shape<0>(xi);
				// evaluate integrand
				auto F = G * N * jac;
				// evaluate singular part
				auto F0 = -N0 * jac0 / (2. * M_PI) * std::log(std::abs(xi(0) - xi0(0)) * jac0);
				
				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [xi0 b];
				xi_t xi = it->get_xi() * ((b(0)-xi0(0))/2.) + (xi0+b)/2.;
				double w = it->get_w() * ((b(0)-xi0(0))/2.);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y);
				// multiply Shape function
				auto N = trial_shape_t::template eval_shape<0>(xi);
				// evaluate integrand
				auto F = G * N * jac;
				// evaluate integrand's singular part
				auto F0 = -N0 * jac0 / (2. * M_PI) * std::log(std::abs(xi(0) - xi0(0)) * jac0);
				
				result.row(i) += (F - F0) * w;
			}
			
			// add analytic integral of singular part
			double d1 = std::abs(a(0) - xi0(0)) * jac0;
			double d2 = std::abs(b(0) - xi0(0)) * jac0;
			result.row(i) += -N0 / (2. * M_PI) * (
				d1 * (std::log(d1) - 1.) + d2 * (std::log(d2) - 1.)
			);
		}
		
		return result;
	}
};


/** \brief Collocational singular integral of the 2D Helmholtz SLP kernel over a constant line element */
template <unsigned expansion_length>
class helmholtz_2d_SLP_collocation_constant_line
{
public:
    /**
     * \brief Evaluate the integral
     * \tparam wavenumber_t the wave number type
     * \param [in] elem the line element
     * \param [in] x0 the singular point
     * \param [in] k the wave number
     * \return the integral value
     */
	template <class wavenumber_t>
	static std::complex<double> eval(
		line_1_elem const &elem,
		line_1_elem::x_t const &x0,
		wavenumber_t const &k)
	{
		double const eulergamma = 0.57721566490153286060;
		std::complex<double> const c(eulergamma, M_PI / 2.);

		// compute elem radius
		auto R = (x0 - elem.get_coords().col(0)).norm();

		auto Q = k * R / 2.;
		auto clnq = c + std::log(Q);

		auto res(clnq - 1.);	// initial (k=0) value of the result
		decltype(Q) B(1.);		// the power term in the series
		double Cn = 0.0;		// the harmonic series up to n = 0
		for (unsigned n = 1; n <= expansion_length; ++n)
		{
			B *= -Q*Q / n / n;		// the actual power term
			Cn += 1. / n;		// the actual harmonic term
			res += B / (2 * n + 1) * (clnq - Cn - 1. / (2 * n + 1));
		}

		return -R / M_PI * res;
	}
};



/** \brief Galerkin singular integral of the 2D Helmholtz SLP kernel over a constant line element */
template <unsigned expansion_length>
class helmholtz_2d_SLP_galerkin_constant_line
{
public:
	template <class wavenumber_t>
	static std::complex<double> eval(double R, wavenumber_t const &k)
	{
		double const eulergamma = 0.57721566490153286060;
		wavenumber_t kR = k * R; 
		wavenumber_t logkR = std::log(kR);
		
		std::complex<double> I = 1. - std::complex<double>(0., 2./M_PI) * (logkR - 1.5 + eulergamma);
		
		wavenumber_t q = -kR * kR;
		wavenumber_t pow = 1.;
		double Cn = 0.;
		
		for (unsigned n = 1; n <= expansion_length; ++n)
		{
			unsigned d = (n+1) * (2*n+1);
			double Fn = 1./d;
			wavenumber_t Gn = logkR/d - (4*n + 3)/2./(d*d);
			pow *= q/(n*n);
			Cn += 1./n;
			I += (Fn-std::complex<double>(0., 2./M_PI)*(Gn+(eulergamma - Cn)*Fn)) * pow;
		}
		
		return I * (R*R) * std::complex<double>(0., -1.);
	}
};


/** \brief Collocational integral of the 2D DLP kernel over a general curved line with general shape sets */
template <class TestField, class TrialField, size_t order>
class helmholtz_2d_DLP_collocation_general
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	
	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;
	
	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;
	
	typedef Eigen::Matrix<std::complex<double>, nTest, nTrial> result_t;
	
	typedef typename trial_field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;
	
	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	
	typedef regular_quad_store<domain_t, order> quadr_t;
	
public:
	
	template <class wave_number_t>
	static result_t eval(elem_t const &elem, wave_number_t const &k)
	{
		result_t result = result_t::Zero();
		
		helmholtz_2d_DLP_kernel<wave_number_t> kernel(k);
		
		xi_t const &a = domain_t::get_corner(0);
		xi_t const &b = domain_t::get_corner(1);
		
		// traverse collocation points
		for (size_t i = 0; i < nTest; ++i)
		{
			xi_t xi0 = test_shape_t::corner_at(i);
			
			// singular point and normal at singular point
			x_t x = elem.get_x(xi0);
			x_t Jxvec = elem.get_normal(xi0);
			x_t nx = Jxvec.normalized();
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [a xi0]
				xi_t xi = it->get_xi() * ((xi0(0)-a(0))/2.) + (xi0+a)/2.;
				double w = it->get_w() * ((xi0(0)-a(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				
				// integrate difference numerically
				result.row(i) += F * w;
			}
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [xi0 b];
				xi_t xi = it->get_xi() * ((b(0)-xi0(0))/2.) + (xi0+b)/2.;
				double w = it->get_w() * ((b(0)-xi0(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				
				// integrate difference numerically
				result.row(i) += F * w;
			}
		}
		return result;
	}
};


/** \brief Collocational integral of the 2D DLPt kernel over a general curved line with general shape sets */
template <class TestField, class TrialField, size_t order>
class helmholtz_2d_DLPt_collocation_general
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	
	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;
	
	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;
	
	typedef Eigen::Matrix<std::complex<double>, nTest, nTrial> result_t;
	
	typedef typename trial_field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;
	
	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	
	typedef regular_quad_store<domain_t, order> quadr_t;
	
public:
	template <class wave_number_t>
	static result_t eval(elem_t const &elem, wave_number_t const &k)
	{
		result_t result = result_t::Zero();
		
		helmholtz_2d_DLPt_kernel<wave_number_t> kernel(k);
		
		xi_t const &a = domain_t::get_corner(0);
		xi_t const &b = domain_t::get_corner(1);
		
		// traverse collocation points
		for (size_t i = 0; i < nTest; ++i)
		{
			xi_t xi0 = test_shape_t::corner_at(i);
			
			// singular point and normal at singular point
			x_t x = elem.get_x(xi0);
			x_t Jxvec = elem.get_normal(xi0);
			x_t nx = Jxvec.normalized();
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [a xi0]
				xi_t xi = it->get_xi() * ((xi0(0)-a(0))/2.) + (xi0+a)/2.;
				double w = it->get_w() * ((xi0(0)-a(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				
				// integrate difference numerically
				result.row(i) += F * w;
			}
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [xi0 b];
				xi_t xi = it->get_xi() * ((b(0)-xi0(0))/2.) + (xi0+b)/2.;
				double w = it->get_w() * ((b(0)-xi0(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				
				// integrate difference numerically
				result.row(i) += F * w;
			}
		}
		return result;
	}
};


/** \brief Collocational integral of the 2D HSP kernel over a general curved line with general shape sets
 * Full singularity subtraction in the reference coordinate system.
 * The singularpart is integrated analytically in HFP sense.
 * The regular part is integrated numerically with standard Gaussian quadrature.
 */
template <class TestField, class TrialField, size_t order>
class helmholtz_2d_HSP_collocation_general
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	
	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;
	
	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;
	
	typedef Eigen::Matrix<std::complex<double>, nTest, nTrial> result_t;
	
	typedef typename trial_field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;
	
	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	
	typedef regular_quad_store<domain_t, order> quadr_t;
	
public:
	
	/** \brief evaluate the singular integral
	 * \param [in] elem the element
	 * \param [in] k the wave number
	 */
	template <class WaveNumber>
	static result_t eval(elem_t const &elem, WaveNumber  const &k)
	{
		// instantiate and clear result matrix
		result_t result = result_t::Zero();
		
		// instantiate the kernel
		helmholtz_2d_HSP_kernel<WaveNumber> kernel(k);
		
		// integration limits in the reference domain
		xi_t const &a = domain_t::get_corner(0);
		xi_t const &b = domain_t::get_corner(1);
		
		// traverse collocation points
		for (size_t i = 0; i < nTest; ++i)
		{
			// the collocation point
			xi_t const &xi0 = test_shape_t::corner_at(i);
			
			// trial shape function and its derivative at the singular point
			auto N0 = trial_shape_t::template eval_shape<0>(xi0);
			auto N1 = trial_shape_t::template eval_shape<1>(xi0);
			
			// singular point and normal at singular point in physical domain
			x_t x = elem.get_x(xi0);
			x_t Jxvec = elem.get_normal(xi0);
			x_t nx = Jxvec.normalized();
			double jac0 = Jxvec.norm();
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [a xi0]
				xi_t xi = it->get_xi() * ((xi0(0)-a(0))/2.) + (xi0+a)/2.;
				double w = it->get_w() * ((xi0(0)-a(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				// evaluate singular part
				auto F0 = (
					(N0/rho + N1)/rho/jac0
					- k*k/2. * N0 * jac0 * std::log(std::abs(rho) * jac0)
				) / (2. * M_PI);
				
				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [xi0 b];
				xi_t xi = it->get_xi() * ((b(0)-xi0(0))/2.) + (xi0+b)/2.;
				double w = it->get_w() * ((b(0)-xi0(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				// evaluate singular part
				auto F0 = (
					(N0/rho + N1)/rho/jac0
					- k*k/2. * N0 * jac0 * std::log(std::abs(rho) * jac0)
				) / (2. * M_PI);
				
				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}
			
			// add analytic integral of singular part
			double d1 = std::abs(a(0) - xi0(0)) * jac0;
			double d2 = std::abs(b(0) - xi0(0)) * jac0;
			result.row(i) += (
				((1./(a(0)-xi0(0)) - 1./(b(0)-xi0(0))) * N0 +
				std::log(std::abs((b(0)-xi0(0)) / (a(0)-xi0(0)))) * N1) / jac0
				- k*k/2. *
				N0 * (d1 * (std::log(d1) - 1.) + d2 * (std::log(d2) - 1.))
				) / (2. * M_PI);
		}
		return result;
	}
};



/** \brief Collocational singular integral of the 2D Helmholtz HSP kernel over a straight line element */
template <class TestField, class TrialField, size_t order>
class helmholtz_2d_HSP_collocation_straight_line
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	
	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;
	
	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;
	
	typedef Eigen::Matrix<std::complex<double>, nTest, nTrial> result_t;
	
	typedef NiHu::line_1_elem elem_t;
	typedef typename elem_t::domain_t domain_t;
	
	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	
	typedef regular_quad_store<domain_t, order> quadr_t;
	
public:
	
	/** \brief evaluate the singular integral
	 * \param [in] elem the element
	 * \param [in] k the wave number
	 */
	template <class WaveNumber>
	static result_t eval(elem_t const &elem, WaveNumber  const &k)
	{
		// instantiate and clear result matrix
		result_t result = result_t::Zero();
		
		// instantiate the kernel
		helmholtz_2d_HSP_kernel<WaveNumber> kernel(k);
		
		// integration limits in the reference domain
		xi_t const &a = domain_t::get_corner(0);
		xi_t const &b = domain_t::get_corner(1);
		
		// traverse collocation points
		for (size_t i = 0; i < nTest; ++i)
		{
			// the collocation point
			xi_t const &xi0 = test_shape_t::corner_at(i);
			
			// trial shape function and its derivative at the singular point
			auto N0 = trial_shape_t::template eval_shape<0>(xi0);
			auto N1 = trial_shape_t::template eval_shape<1>(xi0);
			
			// singular point and normal in physical domain
			x_t x = elem.get_x(xi0);
			x_t Jvec = elem.get_normal(xi0);
			x_t n = Jvec.normalized();
			double jac = Jvec.norm();
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [a xi0]
				xi_t xi = it->get_xi() * ((xi0(0)-a(0))/2.) + (xi0+a)/2.;
				double w = it->get_w() * ((xi0(0)-a(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, n, n);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				// evaluate singular part
				auto F0 = (
					(N0/rho + N1)/rho/jac
					- k*k/2. * N0 * jac * std::log(std::abs(rho) * jac)
				) / (2. * M_PI);
				
				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}
			
			// traverse quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// transform quadrature to [xi0 b];
				xi_t xi = it->get_xi() * ((b(0)-xi0(0))/2.) + (xi0+b)/2.;
				double w = it->get_w() * ((b(0)-xi0(0))/2.);
				double rho = xi(0) - xi0(0);
				
				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				
				// evaluate Green's function
				std::complex<double> G = kernel(x, y, n, n);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = (G * jac) * N;
				// evaluate singular part
				auto F0 = (
					(N0/rho + N1)/rho/jac
					- k*k/2. * N0 * jac * std::log(std::abs(rho) * jac)
				) / (2. * M_PI);
				
				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}
			
			// add analytic integral of singular part
			double d1 = std::abs(a(0) - xi0(0)) * jac;
			double d2 = std::abs(b(0) - xi0(0)) * jac;
			result.row(i) += (
				((1./(a(0)-xi0(0)) - 1./(b(0)-xi0(0))) * N0 +
				std::log(std::abs((b(0)-xi0(0)) / (a(0)-xi0(0)))) * N1) / jac
				- k*k/2. *
				N0 * (d1 * (std::log(d1) - 1.) + d2 * (std::log(d2) - 1.))
				) / (2. * M_PI);
		}
		return result;
	}
};



/** \brief Collocational singular integral of the 3D Helmholtz SLP kernel over a constant planar element
 * \tparam order the quadrature order of the regular (dynamic) part
 */
template <unsigned order>
class helmholtz_3d_SLP_collocation_constant_plane
{
private:
	/** \brief the Dynamic part of the SLP Green's function
	 * This function returns the regular (dynamic) part of the Green's function
	 * defined as exp(-ikr)/r - 1/r
	 * note that the factor 1/(4pi) is not computed for simplification
	 */
	template <class wavenumber_t>
	static std::complex<double> dynamic_part(double const &r, wavenumber_t const &k)
	{
		std::complex<double> const I(0.0, 1.0);
		return -I*k * std::exp(-I*k*r / 2.0) * sinc(k*r / 2.0);
	}

public:
    /**
     * \brief Evaluate the integral
     * \tparam wavenumber_t the wave number type
     * \param [in] elem the line element
     * \param [in] x0 the singular point ( must be internal ! )
     * \param [in] k the wave number
     * \return the integral value
     */
	template <class elem_t, class wavenumber_t>
	static std::complex<double> eval(
		elem_t const &elem,
		typename elem_t::x_t const &x0,
		wavenumber_t const &k)
	{
		typedef regular_quad_store<typename elem_t::domain_t, order> quadr_t;

		enum { N = elem_t::domain_t::num_corners };

		double r[N], theta[N], alpha[N];
		plane_element_helper(elem, x0, r, theta, alpha);

		// integrate static part analytically
		double I_stat = 0.0;
		for (unsigned i = 0; i < N; ++i)
			I_stat += r[i] * std::sin(alpha[i]) *
			std::log(std::tan((alpha[i] + theta[i]) / 2.0) / std::tan(alpha[i] / 2.0)
			);

		// integrate dynamic part
		std::complex<double> I_dyn = 0.0;
		for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
		{
			double r = (elem.get_x(it->get_xi()) - x0).norm();
			double jac = elem.get_normal(it->get_xi()).norm();
			I_dyn += dynamic_part(r, k) * it->get_w() * jac;
		}

		// assemble result from static and dynamic parts
		return (I_stat + I_dyn) / (4.0 * M_PI);
	}
};

/** \brief Collocational singular integral of the 3D Helmholtz HSP kernel over a constant planar element */
template <unsigned order>
class helmholtz_3d_HSP_collocation_constant_plane
{
private:
	template <class wavenumber_t>
	static std::complex<double> dynamic_part(double const &r, wavenumber_t const &k)
	{
		std::complex<double> const I(0.0, 1.0);	// imaginary unit
		std::complex<double> const ikr(I*k*r);
		if (std::abs(r) > 1e-3)
			return (std::exp(-ikr)*(1.0 + ikr) - 1.0 + ikr*ikr / 2.0) / r / r / r;
		else
			return -I*k*k*k * (
			1.0 / 3.0 - ikr*(1.0 / 8.0 - ikr*(1.0 / 30.0 - ikr*(1.0 / 144.0 - ikr*(1.0 / 840.0 - ikr / 5760.0))))
			);
	}

public:
    /**
     * \brief Evaluate the integral
     * \tparam wavenumber_t the wave number type
     * \param [in] elem the line element
     * \param [in] x0 the singular point
     * \param [in] k the wave number
     * \return the integral value
     */
	template <class elem_t, class wavenumber_t>
	static std::complex<double> eval(
		elem_t const &elem,
		typename elem_t::x_t const &x0,
		wavenumber_t const &k)
	{
		
		typedef regular_quad_store<typename elem_t::domain_t, order> quadr_t;
		enum { N = elem_t::domain_t::num_corners };

#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("Integrating Helmholtz HSP over constant plane, N = %d\n", N);
			printed = true;
		}
#endif

		double r[N], theta[N], alpha[N];
		plane_element_helper(elem, x0, r, theta, alpha);

		// integrate static part
		double IG0 = 0.0, IddG0 = 0.0;
		for (unsigned i = 0; i < N; ++i)
		{
			// integral is zero for highly distorted elements
			/** \todo replace with Taylor series expansion */
			if (theta[i] < 1e-3)
				continue;
			IG0 += r[i] * std::sin(alpha[i]) *
				std::log(std::tan((alpha[i] + theta[i]) / 2.0) / tan(alpha[i] / 2.0));
			IddG0 += (std::cos(alpha[i] + theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));
		}

		// integrate dynamic part
		std::complex<double> I_acc = 0.0;
		for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
		{
			double r = (elem.get_x(it->get_xi()) - x0).norm();
			double jac = elem.get_normal(it->get_xi()).norm();
			I_acc += dynamic_part(r, k) * it->get_w() * jac;
		}
		
		// assemble result from static and dynamic parts
		return (IddG0 + k*k / 2.0 * IG0 + I_acc) / (4.0 * M_PI);
	}
};

/** \brief Collocational singular integral of the 2d SLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_SLP_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_SLP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0, 0) = helmholtz_2d_SLP_collocation_constant_line<9>::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center(),
			kernel.derived().get_wave_number());

		return result;
	}
};


/** \brief Collocational singular integral of the 2d SLP kernel over not a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_SLP_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	!(std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value)
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_SLP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result = helmholtz_2d_SLP_collocation_general<TestField, TrialField, 30>::eval(
			trial_field.get_elem(),
			kernel.derived().get_wave_number());

		return result;
	}
};


/** \brief Galerkin singular integral of the 2d SLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_SLP_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value &&
	std::is_same<typename TestField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_SLP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		auto const &elem = trial_field.get_elem();
		double R = (elem.get_coords().col(1) - elem.get_coords().col(0)).norm()/2.;
		result(0, 0) = helmholtz_2d_SLP_galerkin_constant_line<5>::eval(
			R, kernel.derived().get_wave_number());

		return result;
	}
};



/** \brief Collocational singular integral of the 2D DLP kernel over a curved line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_DLP_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	!(std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value)
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_DLP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result = helmholtz_2d_DLP_collocation_general<TestField, TrialField, 30>::eval(
			trial_field.get_elem(),
			kernel.derived().get_wave_number());
			
		return result;
	}
};


/** \brief Collocational singular integral of the 2D DLPt kernel over a curved line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_DLPt_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	!(std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value)
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_DLPt_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result = helmholtz_2d_DLPt_collocation_general<TestField, TrialField, 30>::eval(
			trial_field.get_elem(),
			kernel.derived().get_wave_number());
			
		return result;
	}
};




/** \brief Collocational singular integral of the 2D HSP kernel over a straight line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_HSP_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_HSP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result = helmholtz_2d_HSP_collocation_straight_line<TestField, TrialField, 30>::eval(
			trial_field.get_elem(),
			kernel.derived().get_wave_number());
			
		return result;
	}
};


/** \brief Collocational singular integral of the 2D HSP kernel over a curved line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_HSP_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	!std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_HSP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result = helmholtz_2d_HSP_collocation_general<TestField, TrialField, 30>::eval(
			trial_field.get_elem(),
			kernel.derived().get_wave_number());
			
		return result;
	}
};



/** \brief Collocational singular integral of the 3D SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_SLP_kernel<WaveNumber>, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_3d_SLP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0, 0) = helmholtz_3d_SLP_collocation_constant_plane<7>::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center(),
			kernel.derived().get_wave_number());
		return result;
	}
};


/** \brief Collocational singular integral of the 3D HSP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_HSP_kernel<WaveNumber>, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_3d_HSP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0, 0) = helmholtz_3d_HSP_collocation_constant_plane<7>::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center(),
			kernel.derived().get_wave_number());
		return result;
	}
};

/** \brief Collocational singular integral of the 3D HSP kernel NOT over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_HSP_kernel<WaveNumber>, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		!(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] trial_field the trial and test field
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_3d_HSP_kernel<WaveNumber> > const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("Calling Guiggiani now\n");
			printed = true;
		}
#endif
		typedef guiggiani<TrialField, helmholtz_3d_HSP_kernel<WaveNumber>, 5, 9> guiggiani_t;
		auto const &elem = trial_field.get_elem();
		guiggiani_t gui(elem, kernel.derived());
		for (unsigned r = 0; r < TestField::num_dofs; ++r)
		{
			auto const &xi0 = TestField::nset_t::corner_at(r);
			gui.integrate(result.row(r), xi0, elem.get_normal(xi0));
		}
		
		return result;
	}
};

} // end of namespace NiHu


#endif // NIHU_HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED

