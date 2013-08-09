/** \file poisson_singular_integrals.hpp
 * \brief Analytical expressions for the singular integrals of Poisson kernels
 * \details analytical expression for the Poisson kernels over plane triangles
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED
#define POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../bem/integral_operator.hpp"
#include "poisson_kernel.hpp"

/** \brief Collocational singular integrals of the DLP and DLPt kernels over a plane elements
 * \tparam Formalism the integration formalism (collocational or general)
 * \tparam Kernel the kernel type the test field type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class Kernel, class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational, Kernel, TestField, TrialField,
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
		result_t &result, Kernel const &, TestField const &, TrialField const &)
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
	formalism::collocational, poisson_2d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, line_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		poisson_2d_SLP_kernel const &,
		TestField const &test_field,
		TrialField const &trial_field)
	{
		auto const &c = test_field.get_elem().get_center();
		auto const &C = trial_field.get_elem().get_coords();
		double d1 = (c - C.col(0)).norm();
		double d2 = (c - C.col(1)).norm();
		double d = (C.col(1) - C.col(0)).norm();
		result(0,0) = (d - d1*std::log(d1) - d2*std::log(d2)) / (2*M_PI);
		return result;
	}
};


/** helper function to compute angles and radius in a plane triangle */
template <class scalar_t>
void planar_triangle_helper(
	tria_1_elem const &elem,
	scalar_t r[],
	scalar_t theta[],
	scalar_t alpha [])
{
	auto const &C_old = elem.get_coords();
	auto const &x0 = elem.get_center();
	unsigned const N = tria_1_elem::num_nodes;

	typename tria_1_elem::coords_t R, C;
	for (unsigned i = 0; i < N; ++i)
	{
		R.col(i) = C_old.col(i) - x0;
		r[i] = R.col(i).norm();
		R.col(i) /= r[i];
		C.col(i) = C_old.col(i) - C_old.col((i+1) % N);
		C.col(i) /= C.col(i).norm();
	}

	for (unsigned i = 0; i < N; ++i)
	{
		theta[i] = std::acos(R.col(i).dot(R.col((i+1) % N)));
		alpha[i] = std::acos(R.col(i).dot(C.col(i)));
	}
}


/** \brief collocational singular integral of the SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational, poisson_3d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		poisson_3d_SLP_kernel const &,
		TestField const &,
		TrialField const &trial_field)
	{
		unsigned const N = tria_1_elem::num_nodes;
		double r[N], theta[N], alpha[N];
		planar_triangle_helper(trial_field.get_elem(), r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result(0,0) += r[i] * std::sin(alpha[i]) * std::log(std::tan((alpha[i]+theta[i])/2.0)/tan(alpha[i]/2.0));

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
	formalism::collocational, poisson_3d_HSP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		poisson_3d_HSP_kernel const &,
		TestField const &,
		TrialField const &trial_field)
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
