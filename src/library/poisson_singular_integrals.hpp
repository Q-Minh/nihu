/** \file poisson_singular_integrals.hpp
 * \brief Analytical expressions for the singular integrals of Poisson kernels
 * \details analytical expression for the Poisson kernels over plane triangles
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED
#define POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../bem/integral_operator.hpp"
#include "poisson_kernel.hpp"
#include "plane_triangle_helper.hpp"

/** \brief Singular integrals of the DLP and DLPt kernels over a plane elements
 * \tparam Formalism the integration formalism (collocational or general)
 * \tparam Kernel the kernel type the test field type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class Kernel, class TestField, class TrialField>
class singular_integral_shortcut<
	Kernel, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		(
		( std::is_same<Kernel, poisson_2d_DLP_kernel>::value ||
		  std::is_same<Kernel, poisson_2d_DLPt_kernel>::value
		) && std::is_same<typename TrialField::lset_t, line_1_shape_set>::value
		) || (
		( std::is_same<Kernel, poisson_3d_DLP_kernel>::value ||
		  std::is_same<Kernel, poisson_3d_DLPt_kernel>::value
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
	constexpr static result_t &eval(
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
	poisson_2d_SLP_kernel, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<poisson_2d_SLP_kernel> const &,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		auto const &c = test_field.get_elem().get_center();
		auto const &C = trial_field.get_elem().get_coords();
		auto d1 = (c - C.col(0)).norm();
		auto d2 = (c - C.col(1)).norm();
		result(0,0) = (d1 * (1.0 - std::log(d1)) + d2 * (1.0 - std::log(d2))) / (2*M_PI);
		return result;
	}
};


/** \brief collocational singular integral of the 2D HSP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	poisson_2d_HSP_kernel, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<poisson_2d_HSP_kernel> const &,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		auto const &c = test_field.get_elem().get_center();
		auto const &C = trial_field.get_elem().get_coords();
		auto d1 = (c - C.col(0)).norm();
		auto d2 = (c - C.col(1)).norm();
		result(0,0) = -(1.0/d1 + 1.0/d2) / (2.0*M_PI);
		return result;
	}
};


/** \brief collocational singular integral of the SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	poisson_3d_SLP_kernel, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<poisson_3d_SLP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		unsigned const N = tria_1_elem::num_nodes;
		double r[N], theta[N], alpha[N];
		planar_triangle_helper(trial_field.get_elem(), r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result(0,0) += r[i] * std::sin(alpha[i]) * std::log(std::tan((alpha[i]+theta[i])/2.0)/std::tan(alpha[i]/2.0));

		result(0,0) /= (4.0 * M_PI);

		return result;
	}
};


/** \brief collocational singular integral of the HSP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	poisson_3d_HSP_kernel, TestField, TrialField, singularity::face_match_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<poisson_3d_HSP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		unsigned const N = tria_1_elem::num_nodes;

		double r[N], theta[N], alpha[N];	// radius lengths
		planar_triangle_helper(trial_field.get_elem(), r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result(0,0) += (std::cos(alpha[i]+theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));

		result(0,0) /= (4.0 * M_PI);

		return result;
	}
};

#endif // POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED
