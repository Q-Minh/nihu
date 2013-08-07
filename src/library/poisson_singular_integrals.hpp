/** \file poisson_singular_integrals.hpp
 * \brief Analytical expressions for the singular integrals of Poisson kernels
 * \details analytical expression for the Poisson kernels over plane triangles
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED
#define POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../bem/integral_operator.hpp"
#include "poisson_kernel.hpp"

/** \brief collocational singular integral of the SLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational,
	poisson_2d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		poisson_3d_SLP_kernel const &,
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


/** \brief collocational singular integral of the DLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
/*
template <class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational,
	poisson_2d_DLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<
			typename TrialField::nset_t,
			line_0_shape_set
		>::value
	>::type
>
{
public:
	template <class result_t>
	constexpr static result_t &eval(
		result_t &result,
		poisson_2d_DLP_kernel const &,
		TestField const &test_field,
		TrialField const &trial_field)
	{
		return result;
	}
};
*/


/** \brief collocational singular integral of the SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational,
	poisson_3d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		poisson_3d_SLP_kernel const &,
		TestField const &test_field,
		TrialField const &trial_field)
	{
		auto const &C_old = trial_field.get_elem().get_coords();
		auto const &x0 = test_field.get_elem().get_center();
		unsigned const N = tria_1_elem::num_nodes;

		double r[N];	// radius lengths
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
			double theta = std::acos(R.col(i).dot(R.col((i+1) % N)));
			double alpha = std::acos(R.col(i).dot(C.col(i)));
			result(0,0) += r[i] * std::sin(alpha) * std::log(std::tan((alpha+theta)/2.0)/tan(alpha/2.0));
		}

		result(0,0) /= (4.0 * M_PI);

		return result;
	}
};


/** \brief collocational singular integral of the DLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational,
	poisson_3d_DLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	constexpr static result_t &eval(
		result_t &result,
		poisson_3d_DLP_kernel const &,
		TestField const &,
		TrialField const &)
	{
		return result;
	}
};

#endif // POISSON_SINGULAR_INTEGRALS_HPP_INCLUDED
