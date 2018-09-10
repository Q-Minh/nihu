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

/** \file matsumoto_2010.hpp
 * \brief Explicite hypersingular integrals for collocational Helmholtz BEM with constant triangles
 * \details This file contains the implementation of the collocational integrals for the
 * Burton-Miller formulation over constant triangular elements
 * \ingroup library
 */
#ifndef MATSUMOTO_2010_HPP_INCLUDED
#define MATSUMOTO_2010_HPP_INCLUDED

#include "../core/element_match.hpp"
#include "../core/integral_operator.hpp"
#include "../core/singular_integral_shortcut.hpp"
#include "../core/match_types.hpp"
#include "helmholtz_kernel.hpp"
#include "lib_element.hpp"


namespace NiHu
{

/** \brief internal namespace hiding the stored line quadrature */
namespace matsumoto_internal
{

	/** \brief store-wrapper of a statically stored line quadrature */
	template <unsigned order>
	struct line_quad_store
	{
		/** \brief the stored static quadrature member */
		static gaussian_quadrature<line_domain> const quadrature;
	};

	/** \brief definition of the statically stored quadrature member */
	template <unsigned order>
	gaussian_quadrature<line_domain> const line_quad_store<order>::quadrature(order);

}

/** \brief Collocational hypersingular integral of the HSP kernel over a constant triangle
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
	/** \brief evaluate the singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result the result matrix reference
	 * \param [in] kernel reference to the kernel to be integrated
	 * \param [in] trial_field the trial field instance
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_3d_HSP_kernel<WaveNumber> > const & kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		// Define the imaginary unit
		std::complex<scalar_t> const I(0.0, 1.0);
		// Get wave number from the kernel
		auto const &k = kernel.derived().get_wave_number();

		// Obtain the element
		auto const &tr_elem = trial_field.get_elem();
		unsigned const N = tria_1_elem::num_nodes;

		// Obtain the coordinates and the center
		auto const &coords = tr_elem.get_coords();
		auto const &x0 = tr_elem.get_center();

		for(unsigned i = 0; i < N; ++i)
		{
			// Obtain the nodal distances
			x_t const D = coords.col((i+1)%N) - coords.col(i);
			double l = D.norm();

			// Iterate through quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// Actual point
				x_t const x1 = coords.col((i+1)%N) * (1.0 + it->get_xi()(0))/2
							 + coords.col(i) * (1.0 - it->get_xi()(0))/2;
				// Calculate distance from quadrature point
				x_t const R = x1 - x0;
				double r = R.norm();
				// Calculate the Jacobian
				double tmp = R.dot(D) / (r*l);
				double jac = sqrt(1.0 - tmp*tmp) / r * l / 2.0;
				// Accumulate the result
				result(0,0) += jac * it->get_w() * (std::exp(-I*k*r) / r);
			}
		}
		result(0,0) /= (-4 * M_PI);
		result(0,0) -= I*k / 2.0;
		return result;
	}
private:
	enum { quadrature_order = 31 };
	typedef matsumoto_internal::line_quad_store<quadrature_order> quadr_t;
	typedef typename tria_1_elem::xi_t xi_t;
	typedef typename tria_1_elem::x_t x_t;
	typedef typename kernel_base<helmholtz_3d_HSP_kernel<WaveNumber> >::scalar_t scalar_t;
};


/** \brief Collocational singular integral of the SLP kernel over a constant triangle
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
	/** \brief evaluate the singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result the result matrix reference
	 * \param [in] kernel the kernel to be integrated
	 * \param [in] trial_field the trial field instance
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_3d_SLP_kernel<WaveNumber> > const & kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		// Define the imaginary unit
		std::complex<scalar_t> const I(0.0, 1.0);
		// Get wavenumber from the kernel
		auto const &k = kernel.derived().get_wave_number();

		// Obtain the element
		auto const &tr_elem = trial_field.get_elem();
		unsigned const N = tria_1_elem::num_nodes;

		// Obtain the coordinates and the center
		auto const &coords = tr_elem.get_coords();
		auto const &x0 = tr_elem.get_center();

		for(unsigned i = 0; i < N; ++i)
		{
			// Obtain the nodal distances
			x_t const D = coords.col((i+1)%N) - coords.col(i);
			double l = D.norm();


			// Iterate through quadrature points
			for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			{
				// Actual point
				x_t const x1 = coords.col((i+1)%N) * (1.0 + it->get_xi()(0))/2
							 + coords.col(i) * (1.0 - it->get_xi()(0))/2;

				// Calculate distance from quadrature point
				x_t const R = x1 - x0;
				double r = R.norm();
				// Calculate the Jacobian
				double tmp = R.dot(D) / (r*l);
				double jac = sqrt(1.0 - tmp*tmp) / r * l / 2.0;
				// Accumulate the result
				result(0,0) += jac * it->get_w() * (std::exp(-I*k*r));
			}
		}
		result(0,0) *= I / (4.0 * M_PI * k);
		result(0,0) -= I/(2*k);
		return result;
	}
private:
	enum { quadrature_order = 31 };
	typedef matsumoto_internal::line_quad_store<quadrature_order> quadr_t;
	typedef typename tria_1_elem::xi_t xi_t;
	typedef typename tria_1_elem::x_t x_t;
	typedef typename kernel_base<helmholtz_3d_SLP_kernel<WaveNumber> >::scalar_t scalar_t;
};

}

#endif // MATSUMOTO_2010_HPP_INCLUDED

