/** \file helmholtz_singular_integrals.hpp
 * \brief (Semi)analytical expressions for the singular integrals of Helmholtz kernels
 * \details Semianalytical expression for the helmholtz kernels over plane triangles
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
#define HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../bem/integral_operator.hpp"
#include "helmholtz_kernel.hpp"

#include "../util/math_functions.hpp"

/** \brief Trivial collocational singular integrals of various kernels over plane surfaces
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class Kernel, class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational,
	Kernel, TestField, TrialField,
	typename std::enable_if<
		( std::is_same<Kernel, helmholtz_3d_DLP_kernel>::value ||
		  std::is_same<Kernel, helmholtz_3d_DLPt_kernel>::value
		) && std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value
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
	static constexpr result_t &eval(
		result_t &result, Kernel const &, TestField const &, TrialField const &)
	{
		return result;
	}
};

/** \brief Collocational singular integral of the SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational, helmholtz_3d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
private:
	/** \brief Compute the regular dynamic part of the singular kernel
	 * \tparam T the scalar type
	 * \param [in] r the scalar distance
	 * \param [in] k the wave number
	 * \return the dynamic part of the singular kernel
	 */
	template <class T>
	static std::complex<T> dynamic_part(T const &r, std::complex<T> const &k)
	{
		std::complex<T> const I(0.0, 1.0);	// imaginary unit
		return -I*k * std::exp(-I*k*r/2.0) * sinc(k*r/2.0);
	}

public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		helmholtz_3d_SLP_kernel const &kernel,
		TestField const &test_field,
		TrialField const &trial_field)
	{
		unsigned const N = tria_1_elem::num_nodes;

		auto const &tr_elem = trial_field.get_elem();
		auto const &C_old = tr_elem.get_coords();
		auto const &x0 = test_field.get_elem().get_center();

		// integrate static part analytically
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

		// integrate dynamic_part
		std::complex<double> I_acc = 0.0;
		for (auto it = m_quadrature.begin(); it != m_quadrature.end(); ++it)
			I_acc += dynamic_part(
				(tr_elem.get_x(it->get_xi()) - x0).norm(),
				kernel.get_wave_number()
			) * it->get_w();
		// multiply by jacobian
		I_acc *= tr_elem.get_normal(gauss_tria::xi_t()).norm();

		result(0,0) += I_acc;

		result(0,0) /= (4.0 * M_PI);

		return result;
	}

private:
	/** \brief Regular quadrature used to integrate the dynamic part */
	static gauss_tria const m_quadrature;
};

/** \brief Static regular quadrature instance */
template <class TestField, class TrialField>
gauss_tria const
singular_integral_shortcut<
	formalism::collocational, helmholtz_3d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>::m_quadrature(7);

#endif // HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
