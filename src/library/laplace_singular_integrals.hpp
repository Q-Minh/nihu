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

/** \file laplace_singular_integrals.hpp
 * \brief Analytical expressions for the singular integrals of laplace kernels
 * \details analytical expression for the laplace kernels over plane triangles
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED
#define LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../core/integral_operator.hpp"
#include "laplace_kernel.hpp"
#include "plane_triangle_helper.hpp"


class laplace_2d_SLP_collocation_constant_line
{
public:
	static double eval(line_1_elem const &elem, line_1_elem::x_t const &x0)
	{
		auto const &C = elem.get_coords();
		double d1 = (x0 - C.col(0)).norm(), d2 = (x0 - C.col(1)).norm();
		return (d1 * (1.0 - std::log(d1)) + d2 * (1.0 - std::log(d2))) / (2.0 * M_PI);
	}
};

class laplace_2d_HSP_collocation_constant_line
{
public:
	static double eval(line_1_elem const &elem, line_1_elem::x_t const &x0)
	{
		auto const &C = elem.get_coords();
		double d1 = (x0 - C.col(0)).norm(), d2 = (x0 - C.col(1)).norm();
		return -(1.0 / d1 + 1.0 / d2) / (2.0 * M_PI);
	}
};


class laplace_3d_SLP_collocation_constant_triangle
{
public:
	static double eval(tria_1_elem const &elem, tria_1_elem::x_t const &x0)
	{
		double r[3], theta[3], alpha[3], result = 0.0;
		planar_triangle_helper(elem, x0, r, theta, alpha);

		for (unsigned i = 0; i < 3; ++i)
			result += r[i] * std::sin(alpha[i]) * std::log(std::tan((alpha[i] + theta[i]) / 2.0) / std::tan(alpha[i] / 2.0));

		return result / (4.0 * M_PI);
	}
};

class laplace_3d_HSP_collocation_constant_triangle
{
public:
	static double eval(tria_1_elem const &elem, tria_1_elem::x_t const &x0)
	{
		double r[3], theta[3], alpha[3], result = 0.0;
		planar_triangle_helper(elem, x0, r, theta, alpha);

		for (unsigned i = 0; i < 3; ++i)
			result += (std::cos(alpha[i] + theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));

		return result / (4.0 * M_PI);
	}
};




/** \brief Singular integrals of the DLP and DLPt kernels over plane triangle elements
 * \tparam Formalism the integration formalism (collocational or general)
 * \tparam Kernel the kernel type the test field type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class Kernel, class TestField, class TrialField>
class singular_integral_shortcut<
	Kernel, TestField, TrialField, match::face_match_type,
	typename std::enable_if<
		(
		( std::is_same<Kernel, laplace_2d_DLP_kernel>::value ||
		  std::is_same<Kernel, laplace_2d_DLPt_kernel>::value
		) && std::is_same<typename TrialField::lset_t, line_1_shape_set>::value
		) || (
		( std::is_same<Kernel, laplace_3d_DLP_kernel>::value ||
		  std::is_same<Kernel, laplace_3d_DLPt_kernel>::value
		) && std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value
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
	CONSTEXPR static result_t &eval(
		result_t &result,
		kernel_base<Kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};


/** \brief collocational singular integral of the 2D SLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::face_match_type,
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
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_SLP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0,0) = laplace_2d_SLP_collocation_constant_line::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center());
		return result;
	}
};


/** \brief Galerkin face-match singular integral of the 2D SLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::face_match_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_SLP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		auto const &C = trial_field.get_elem().get_coords();
		auto d = (C.col(1) - C.col(0)).norm();	// element length
		result(0, 0) = -d*d*(std::log(d)-1.5) / (2.0*M_PI);
		return result;
	}
};


/** \brief Galerkin face-match singular integral of the 2D SLP kernel over a linear line
* \tparam TestField the test field type
* \tparam TrialField the trial field type
*/
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::face_match_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_1_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	* \tparam result_t the result matrix type
	* \param [in, out] result reference to the result
	* \param [in] trial_field the test and trial fields
	* \return reference to the result matrix
	*/
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_SLP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		auto const &C = trial_field.get_elem().get_coords();
		auto d = (C.col(1) - C.col(0)).norm();	// element length
		auto lnd = std::log(d);
		result <<
			lnd - 1.75, lnd - 1.25,
			lnd - 1.25, lnd - 1.45;
		result *= -d*d / (8.0*M_PI);
		return result;
	}
};

/** \brief collocational singular integral of the 2D HSP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_HSP_kernel, TestField, TrialField, match::face_match_type,
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
	 * \param [in] trial_field the test and trial fields
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_HSP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0,0) = laplace_2d_HSP_collocation_constant_line::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center());
		return result;
	}
};


/** \brief Galerkin face-match singular integral of the 2D HSP kernel over a constant line
* \tparam TestField the test field type
* \tparam TrialField the trial field type
*/
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_HSP_kernel, TestField, TrialField, match::face_match_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	// the condition always evaluates to false
	static_assert(!std::is_same<typename TestField::nset_t, typename TrialField::nset_t>::value,
		"\n\nThe 2D HSP kernel of the Laplace equation can not be integrated over a constant line element in a Galerkin sense.\n\n");
};

/** \brief collocational singular integral of the 3D SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_3d_SLP_kernel, TestField, TrialField, match::face_match_type,
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
	 * \param [in] trial_field the trial and test field
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_SLP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0, 0) = laplace_3d_SLP_collocation_constant_triangle::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center());
		return result;
	}
};


/** \brief collocational singular integral of the HSP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_3d_HSP_kernel, TestField, TrialField, match::face_match_type,
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
	 * \param [in] trial_field the trial and test field
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_HSP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0, 0) = laplace_3d_HSP_collocation_constant_triangle::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center());
		return result;
	}
};

#endif // LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED
