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

/** \file elastostatics_singular_integrals.hpp
 * \brief Analytical expressions for the singular integrals of elastostatic kernels
 */
#ifndef ELASTOSTATICS_SINGULAR_INTEGRALS_HPP_INCLUDED
#define ELASTOSTATICS_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../core/integral_operator.hpp"
#include "elastostatics_kernel.hpp"
#include "guiggiani_1992.hpp"


/** \brief collocational singular integral of the T kernel
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	elastostatics_3d_T_kernel, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value
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
		kernel_base<elastostatics_3d_T_kernel> const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		typedef guiggiani<TrialField, elastostatics_3d_T_kernel, 5, 9> guiggiani_t;
		auto const &elem = trial_field.get_elem();
		guiggiani_t gui(elem, kernel.derived());

		auto const &xi0 = TestField::nset_t::corner_at(0);
		gui.integrate(result, xi0, elem.get_normal(xi0));

		return result;
	}
};


#endif // ELASTOSTATICS_SINGULAR_INTEGRALS_HPP_INCLUDED

