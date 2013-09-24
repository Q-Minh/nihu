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

/** \file matsumoto_2010.hpp
 * \brief Explicite hypersingular integrals for collocational Helmholtz BEM with constant triangles
 * \ingroup library
 */
#ifndef MATSUMOTO_2010_HPP_INCLUDED
#define MATSUMOTO_2010_HPP_INCLUDED

#include "../core/integral_operator.hpp"
#include "helmholtz_kernel.hpp"

/** \brief Collocational hypersingular integral of the HSP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	helmholtz_3d_HSP_kernel, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value &&
		std::is_same<typename TestField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate the singular integral
	 * \tparam result_t the result matrix type
	 * \param result the result matrix reference
	 * \param trial_field the trial field instance
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<helmholtz_3d_HSP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		// TODO: fill this
		return result;
	}
};

#endif // MATSUMOTO_2010_HPP_INCLUDED