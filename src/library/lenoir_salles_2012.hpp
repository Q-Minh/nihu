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

/** \file lenoir_salles_2012.hpp
 * \brief Explicite singular Galerkin integrals of Laplace kernels over plane triangles
 * \ingroup library
 */
#ifndef LENOIR_SALLES_2012_HPP_INCLUDED
#define LENOIR_SALLES_2012_HPP_INCLUDED

#include "../core/integral_operator.hpp"
#include "laplace_kernel.hpp"
#include "lib_element.hpp"


/** \brief Galerkin singular integral of the Laplace SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_3d_SLP_kernel, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value &&
		std::is_same<typename TestField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	/** \brief helper function to compute geometrical parameters of a triangle
	 * \tparam scalar_t the scalar type of the parameters
	 * \param [in] elem the triangle element
	 * \param [out] gamma distance of a corner from the opposite side vector
	 * \param [out] splus signed distance of other corner from the projection point
	 * \param [out] sminus signed distance of other corner from the projection point
	 */
	template <class scalar_t>
	static void gamma_splus_sminus(
		tria_1_elem const &elem,
		scalar_t gamma[],
		scalar_t splus[],
		scalar_t sminus[])
	{
		auto const &C = elem.get_coords();
		unsigned const N = tria_1_elem::num_nodes;

		for (unsigned i = 0; i < N; ++i)
		{
			auto a = C.col(i);
			auto b = C.col((i+1)%N);
			auto c = C.col((i+2)%N);
			auto d = c-b;				// opposite side vector
			auto unit_d = d / d.norm();	// unit vector
			auto x = (a-b).dot(unit_d)*unit_d + b;	// projection
			gamma[i] = (x-a).norm();				// distance from projection
			sminus[i] = (x-b).dot(unit_d);			// signed distance
			splus[i] = (x-c).dot(unit_d);			// signed distance
		}
	}

	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param result the result matrix reference
	 * \param trial_field the trial field instance
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_SLP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		unsigned const N = tria_1_elem::num_nodes;
		auto const &elem = trial_field.get_elem();

		// compute geometrical helper values
		double gamma[N], splus[N], sminus[N];
		gamma_splus_sminus(elem, gamma, splus, sminus);

		// element area
		auto S = elem.get_normal().norm()/2.0;

		// the famous integral expression from Lenoir-Salles
		for (unsigned i = 0; i < N; ++i)
			result(0,0) -= gamma[i] * (
				std::asinh(splus[i]/gamma[i]) - std::asinh(sminus[i]/gamma[i])
			);

		result(0,0) *= 2.0/3.0 * S / (4.0*M_PI);

		return result;
	}
};

#endif // LENOIR_SALLES_2012_HPP_INCLUDED
