// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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
 * \details Semianalytical expression for the Helmholtz kernels over plane triangles
 */
#ifndef HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
#define HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../core/integral_operator.hpp"
#include "helmholtz_kernel.hpp"
#include "plane_triangle_helper.hpp"
#include "../util/math_functions.hpp"

/** \brief Trivial integrals of various kernels over plane surfaces
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <template< class WaveNumber> class Kernel, class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	Kernel<WaveNumber>, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		(
		( std::is_same<Kernel<WaveNumber>, helmholtz_3d_DLP_kernel<WaveNumber> >::value ||
		  std::is_same<Kernel<WaveNumber>, helmholtz_3d_DLPt_kernel<WaveNumber> >::value
		) && std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value
		)
		||
		(
		std::is_same<Kernel<WaveNumber>, helmholtz_2d_DLP_kernel<WaveNumber> >::value
		&& std::is_same<typename TrialField::lset_t, line_1_shape_set>::value
		)
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
	static CONSTEXPR result_t &eval(
		result_t &result,
		kernel_base<Kernel<WaveNumber> > const &,
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
	helmholtz_2d_SLP_kernel<WaveNumber>, TestField, TrialField, singularity::face_match_type,
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
        std::complex<double> const c(0.57721566490153286060, M_PI/2);
        unsigned const N = 5;

		// compute elem radius
		auto const &elem = trial_field.get_elem();
		auto R = (elem.get_center() - elem.get_coords().col(0)).norm();

		auto Q = kernel.get_data().get_wave_number() * R / 2.0;
        auto clnq = c + std::log(Q);

        auto res(clnq-1.0);		// initial (k=0) value of the result
        decltype(Q) B(1.0);		// the power term in the series
        double bn = 0.0;		// the harmonic series up to n = 0
        for (unsigned n = 1; n < N; ++n)
        {
			B *= -Q*Q/n/n;		// the actual power term
			bn += 1.0/n;		// the actual harmonic term
			res += B/(2*n+1) * (clnq - bn - 1.0/(2*n+1));
        }

		result(0,0) = -R / M_PI * res;

		return result;
	}
};


/** \brief store-wrapper of a statically stored quadrature */
template <unsigned order>
struct tria_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<tria_domain> const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <unsigned order>
gaussian_quadrature<tria_domain> const tria_quad_store<order>::quadrature(order);


/** \brief Collocational singular integral of the SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_SLP_kernel<WaveNumber>, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
private:
	enum { quadrature_order = 7 };
	typedef tria_quad_store<quadrature_order> quadr_t;

	/** \brief Compute the regular dynamic part of the singular kernel
	 * \tparam T the scalar type
	 * \param [in] r the scalar distance
	 * \param [in] k the wave number
	 * \return the dynamic part of the singular kernel
	 */
	template <class T>
	static std::complex<T> dynamic_part(T const &r, WaveNumber const &k)
	{
		std::complex<T> const I(0.0, 1.0);	// imaginary unit
		return -I*k * std::exp(-I*k*r/2.0) * sinc(k*r/2.0);
	}

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
		auto const &tr_elem = trial_field.get_elem();
		unsigned const N = tria_1_elem::num_nodes;
		double r[N], theta[N], alpha[N];
		planar_triangle_helper(tr_elem, r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result(0,0) += r[i] * std::sin(alpha[i]) *
				std::log(
					std::tan((alpha[i]+theta[i])/2.0)/
					std::tan(alpha[i]/2.0)
				);

		// integrate dynamic_part
		auto const &x0 = tr_elem.get_center();
		std::complex<double> I_dyn = 0.0;
		for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			I_dyn += dynamic_part(
				(tr_elem.get_x(it->get_xi()) - x0).norm(),
				kernel.get_data().get_wave_number()
			) * it->get_w();
		// multiply by Jacobian
		I_dyn *= tr_elem.get_normal(tria_domain::xi_t()).norm();

		result(0,0) += I_dyn;

		result(0,0) /= (4.0 * M_PI);

		return result;
	}
};


/** \brief Collocational singular integral of the HSP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_HSP_kernel<WaveNumber>, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
private:
	enum { quadrature_order = 7 };
	typedef tria_quad_store<quadrature_order> quadr_t;

	/** \brief Compute the regular dynamic part of the singular kernel
	 * \tparam T the scalar type
	 * \param [in] r the scalar distance
	 * \param [in] k the wave number
	 * \return the dynamic part of the singular kernel
	 */
	template <class T>
	static std::complex<T> dynamic_part(T const &r, WaveNumber const &k)
	{
		std::complex<T> const I(0.0, 1.0);	// imaginary unit
		std::complex<T> const ikr(I*k*r);
		if (std::abs(r) > 1e-3)
			return (std::exp(-ikr)*(1.0+ikr)-1.0+ikr*ikr/2.0)/r/r/r;
		else
			return -I*k*k*k * (
				1.0/3.0 - ikr*(1.0/8.0 - ikr*(1.0/30.0 - ikr*(1.0/144.0 - ikr*(1.0/840.0 - ikr/5760.0))))
			);
	}

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
		unsigned const N = tria_1_elem::num_nodes;
		double r[N], theta[N], alpha[N];
		planar_triangle_helper(trial_field.get_elem(), r, theta, alpha);

		double IG0 = 0.0, IddG0 = 0.0;
		for (unsigned i = 0; i < N; ++i)
		{
			IG0 += r[i] * std::sin(alpha[i]) * std::log(std::tan((alpha[i]+theta[i])/2.0)/tan(alpha[i]/2.0));
			IddG0 += (std::cos(alpha[i]+theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));
		}

		// integrate dynamic_part
		auto const &tr_elem = trial_field.get_elem();
		auto const &x0 = tr_elem.get_center();
		std::complex<double> I_acc = 0.0;
		for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			I_acc += dynamic_part(
				(tr_elem.get_x(it->get_xi()) - x0).norm(),
				kernel.get_data().get_wave_number()
			) * it->get_w();
		// multiply by Jacobian
		I_acc *= tr_elem.get_normal(tria_domain::xi_t()).norm();

		// assemble result from static and dynamic parts
		auto k2p2 = kernel.get_data().get_wave_number()*kernel.get_data().get_wave_number()/2.0;
		result(0,0) += (IddG0 + k2p2 * IG0 + I_acc) / (4.0 * M_PI);

		return result;
	}
};

#endif // HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
