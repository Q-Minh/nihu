/** \file helmholtz_singular_integrals.hpp
 * \brief (Semi)analytical expressions for the singular integrals of Helmholtz kernels
 * \details Semianalytical expression for the Helmholtz kernels over plane triangles
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
#define HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../bem/integral_operator.hpp"

#include "helmholtz_kernel.hpp"
#include "../util/math_functions.hpp"
#include "plane_triangle_helper.hpp"

/** \brief Trivial collocational integrals of various kernels over plane surfaces
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, template< class WaveNumber> class Kernel, class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational,
	Kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		( std::is_same<Kernel<WaveNumber>, helmholtz_3d_DLP_kernel<WaveNumber> >::value ||
		  std::is_same<Kernel<WaveNumber>, helmholtz_3d_DLPt_kernel<WaveNumber> >::value
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
		result_t &result, Kernel<WaveNumber> const &, TestField const &, TrialField const &)
	{
		return result;
	}
};

/** \brief store-wrapper of a statically stored quadrature */
template <unsigned order>
struct tria_quad_store
{
	/** \brief the stored static quadrature member */
	static gauss_tria const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <unsigned order>
gauss_tria const tria_quad_store<order>::quadrature(order);


/** \brief Collocational singular integral of the SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational, helmholtz_3d_SLP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
private:
	enum {
		quadrature_order = 7
	};
	typedef tria_quad_store<quadrature_order> quadr_t;

	/** \brief Compute the regular dynamic part of the singular kernel
	 * \tparam T the scalar type
	 * \param [in] r the scalar distance
	 * \param [in] k the wave number
	 * \return the dynamic part of the singular kernel
	 */
	template <class T>
	static std::complex<T> dynamic_part(T const &r, WaveNumber const &k)
	{
		std::complex<T> const I(0.0, 1.0);	// imaginary unit
		return -I*k * std::exp(-I*k*r/2.0) * sinc(k*r/2.0);
	}

public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		helmholtz_3d_SLP_kernel<WaveNumber> const &kernel,
		TestField const &,
		TrialField const &trial_field)
	{
		unsigned const N = tria_1_elem::num_nodes;
		double r[N], theta[N], alpha[N];
		planar_triangle_helper(trial_field.get_elem(), r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result(0,0) += r[i] * std::sin(alpha[i]) *
				std::log(
					std::tan((alpha[i]+theta[i])/2.0)/
					std::tan(alpha[i]/2.0)
				);


		// integrate dynamic_part
		auto const &tr_elem = trial_field.get_elem();
		auto const &x0 = tr_elem.get_center();
		std::complex<double> I_acc = 0.0;
		for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			I_acc += dynamic_part(
				(tr_elem.get_x(it->get_xi()) - x0).norm(),
				kernel.get_wave_number()
			) * it->get_w();
		// multiply by Jacobian
		I_acc *= tr_elem.get_normal(gauss_tria::xi_t()).norm();

		result(0,0) += I_acc;

		result(0,0) /= (4.0 * M_PI);

		return result;
	}
};


/** \brief Collocational singular integral of the HSP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class WaveNumber, class TestField, class TrialField>
class singular_integral_shortcut<
	formalism::collocational, helmholtz_3d_HSP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename TrialField::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
	>::type
>
{
private:
	enum {
		quadrature_order = 7
	};
	typedef tria_quad_store<quadrature_order> quadr_t;

	/** \brief Compute the regular dynamic part of the singular kernel
	 * \tparam T the scalar type
	 * \param [in] r the scalar distance
	 * \param [in] k the wave number
	 * \return the dynamic part of the singular kernel
	 */
	template <class T>
	static std::complex<T> dynamic_part(T const &r, WaveNumber const &k)
	{
		std::complex<T> const I(0.0, 1.0);	// imaginary unit
		std::complex<T> const ikr(I*k*r);
		if (std::abs(r) > 1e-3)
			return (std::exp(-ikr)*(1.0+ikr)-1.0+ikr*ikr/2.0)/r/r/r;
		else
			return -I*k*k*k * (
				1.0/3.0 - ikr*(1.0/8.0 - ikr*(1.0/30.0 - ikr*(1.0/144.0 - ikr*(1.0/840.0 - ikr/5760.0))))
			);
	}

public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		helmholtz_3d_HSP_kernel<WaveNumber> const &kernel,
		TestField const &,
		TrialField const &trial_field)
	{
		unsigned const N = tria_1_elem::num_nodes;
		double r[N], theta[N], alpha[N];
		planar_triangle_helper(trial_field.get_elem(), r, theta, alpha);

		double IG0 = 0.0, IddG0 = 0.0;
		for (unsigned i = 0; i < N; ++i)
		{
			IG0 += r[i] * std::sin(alpha[i]) * std::log(std::tan((alpha[i]+theta[i])/2.0)/tan(alpha[i]/2.0));
			IddG0 += (std::cos(alpha[i]+theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));
		}

		// integrate dynamic_part
		auto const &tr_elem = trial_field.get_elem();
		auto const &x0 = tr_elem.get_center();
		std::complex<double> I_acc = 0.0;
		for (auto it = quadr_t::quadrature.begin(); it != quadr_t::quadrature.end(); ++it)
			I_acc += dynamic_part(
				(tr_elem.get_x(it->get_xi()) - x0).norm(),
				kernel.get_wave_number()
			) * it->get_w();
		// multiply by Jacobian
		I_acc *= tr_elem.get_normal(gauss_tria::xi_t()).norm();

		// assemble result from static and dynamic parts
		auto k2p2 = kernel.get_wave_number()*kernel.get_wave_number()/2.0;
		result(0,0) += (IddG0 + k2p2 * IG0 + I_acc) / (4.0 * M_PI);

		return result;
	}
};

#endif // HELMHOLTZ_SINGULAR_INTEGRALS_HPP_INCLUDED
