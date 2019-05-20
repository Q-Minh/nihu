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
#include "../core/singular_integral_shortcut.hpp"

#include "elastostatics_kernel.hpp"
#include "guiggiani_1992.hpp"
#include "lib_element.hpp"


namespace NiHu
{
	
/** \brief Galerkin face match integral of 2D Elastostatics U kernel over constant line */
class elastostatics_2d_U_galerkin_face_constant_line
{
	typedef Eigen::Matrix<double, 2, 2> result_t;
public:
	/**
	 * \brief Evaluate the integral
	 * \param [in] elem the line element
	 * \return the integral value
	 */
	static result_t eval(line_1_elem const &elem, double nu)
	{
		auto const &C = elem.get_coords();
		auto rvec = (C.col(1) - C.col(0)).eval();
		auto r = rvec.norm(); 	// elem length
		auto gradr = rvec.normalized();
		
		return r*r * (
			-(3. - 4.*nu) * result_t::Identity() * (std::log(r) - 1.5)
			+ gradr * gradr.transpose()
		) / (8.*M_PI * (1-nu));
	}
};


/** \brief Galerkin face-match singular integral of the 2D U kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	elastostatics_2d_U_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value &&
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
		kernel_base<elastostatics_2d_U_kernel> const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result = elastostatics_2d_U_galerkin_face_constant_line::eval(
			trial_field.get_elem(), kernel.derived().get_poisson_ratio());
		return result;
	}
};




/** \brief Galerkin face-match singular integral of the 2D T kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	elastostatics_2d_T_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TestField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value &&
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
		kernel_base<elastostatics_2d_T_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result.setZero();
		return result;
	}
};


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

} // end of namespace NiHu


#endif // ELASTOSTATICS_SINGULAR_INTEGRALS_HPP_INCLUDED

