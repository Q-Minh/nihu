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

/** \file helmholtz_singular_integrals.hpp
 * \brief (Semi)analytical expressions for the singular integrals of Helmholtz kernels over plane elements
 */
#ifndef NIHU_HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
#define NIHU_HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "library/lib_element.hpp"
#include "../core/integral_operator.hpp"
#include "helmholtz_kernel.hpp"
#include "plane_element_helper.hpp"
#include "../util/math_functions.hpp"
#include "guiggiani_1992.hpp"
#include "../core/singular_integral_shortcut.hpp"

namespace NiHu
{

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
		double const EulerMascheroni = 0.57721566490153286060;
		std::complex<double> const c(EulerMascheroni, M_PI / 2.0);

		// compute elem radius
		auto R = (x0 - elem.get_coords().col(0)).norm();

		auto Q = k * R / 2.0;
		auto clnq = c + std::log(Q);

		auto res(clnq - 1.0);		// initial (k=0) value of the result
		decltype(Q) B(1.0);		// the power term in the series
		double bn = 0.0;		// the harmonic series up to n = 0
		for (unsigned n = 1; n < expansion_length; ++n)
		{
			B *= -Q*Q / n / n;		// the actual power term
			bn += 1.0 / n;		// the actual harmonic term
			res += B / (2 * n + 1) * (clnq - bn - 1.0 / (2 * n + 1));
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
		double const gamma = 0.57721566490153286060;
		double D = 2. * R;
		wavenumber_t kR = k * R; 
		wavenumber_t logkR = std::log(kR);
		
		std::complex<double> I = 1. - std::complex<double>(0., 2./M_PI) * (logkR - 1.5 + gamma);
		
		wavenumber_t q = -kR * kR;
		wavenumber_t pow = 1.0;
		double Cn = 0.0;
		
		for (unsigned n = 1; n <= expansion_length; ++n)
		{
			unsigned d = (n+1) * (2*n+1);
			double Fn = 1/d;
			wavenumber_t Gn = logkR/d - (4*n + 3)/2.0/(d*d);
			pow *= q/(n*n);
			Cn += 1./n;
			I += (Fn-std::complex<double>(0., 2./M_PI)*(Gn+(gamma - Cn)*Fn)) * pow;
		}
		
		return I * (D*D) * std::complex<double>(0, -.25);
	}
};




/** \brief store-wrapper of a statically stored quadrature */
template <class domain_t, unsigned order>
struct domain_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<domain_t> const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <class domain_t, unsigned order>
gaussian_quadrature<domain_t> const domain_quad_store<domain_t, order>::quadrature(order);


/** \brief Collocational singular integral of the 3D Helmholtz SLP kernel over a constant planar element */
template <unsigned order>
class helmholtz_3d_SLP_collocation_constant_plane
{
private:
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
		typedef domain_quad_store<typename elem_t::domain_t, order> quadr_t;

		enum { N = elem_t::domain_t::num_corners };

		double r[N], theta[N], alpha[N];
		plane_element_helper(elem, x0, r, theta, alpha);

		// integrate dynamic part
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
			I_dyn += dynamic_part(r, k) * it->get_w() * elem.get_normal(it->get_xi()).norm();
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
		typedef domain_quad_store<typename elem_t::domain_t, order> quadr_t;
		enum { N = elem_t::domain_t::num_corners };

		double r[N], theta[N], alpha[N];
		plane_element_helper(elem, x0, r, theta, alpha);

		// integrate static part
		double IG0 = 0.0, IddG0 = 0.0;
		for (unsigned i = 0; i < N; ++i)
		{
			IG0 += r[i] * std::sin(alpha[i]) * std::log(std::tan((alpha[i] + theta[i]) / 2.0) / tan(alpha[i] / 2.0));
			IddG0 += (std::cos(alpha[i] + theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));
		}

		// integrate dynamic part
		std::complex<double> I_acc = 0.0;
		for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
		{
			double r = (elem.get_x(it->get_xi()) - x0).norm();
			I_acc += dynamic_part(r, k) * it->get_w() * elem.get_normal(it->get_xi()).norm();
		}

		// assemble result from static and dynamic parts
		return (IddG0 + k*k / 2.0 * IG0 + I_acc) / (4.0 * M_PI);
	}
};


/** \brief Trivial integrals of various 3d kernels over plane surfaces
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <template<class WaveNumber> class Kernel, class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	Kernel<WaveNumber>, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		( std::is_same<Kernel<WaveNumber>, helmholtz_3d_DLP_kernel<WaveNumber> >::value ||
		  std::is_same<Kernel<WaveNumber>, helmholtz_3d_DLPt_kernel<WaveNumber> >::value
		) && std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate the kernel (zero)
	 * \tparam result_t the result's type
	 * \param [in] result the result reference
	 * \return the result reference
	 */
	template <class result_t>
	static constexpr result_t &eval(
		result_t &result,
		kernel_base<Kernel<WaveNumber> > const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};

/** \brief Trivial integrals of various 2d kernels over plane surfaces
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_2d_DLP_kernel<WaveNumber>, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, line_1_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate the kernel (zero)
	 * \tparam result_t the result's type
	 * \param [in] result the result reference
	 * \return the result reference
	 */
	template <class result_t>
	static constexpr result_t &eval(
		result_t &result,
		kernel_base<helmholtz_2d_DLP_kernel<WaveNumber> > const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
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
	std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
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
		result(0, 0) = helmholtz_2d_SLP_collocation_constant_line<5>::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center(),
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
	std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value &&
	std::is_same<typename TestField::lset_t, line_1_shape_set>::value &&
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


/** \brief Collocational singular integral of the 3D SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_SLP_kernel<WaveNumber>, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
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
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
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

/** \brief Collocational singular integral of the 3D HSP kernel not over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_HSP_kernel<WaveNumber>, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		!(std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
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

}


#endif // NIHU_HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED

