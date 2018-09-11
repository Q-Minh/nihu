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

/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for CPV and HPF collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include "location_normal.hpp"
#include "../core/field.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../core/kernel.hpp"
#include "../core/shapeset.hpp"
#include "../util/math_constants.hpp"
#include "../util/block_product.hpp"
#include "../util/store_pattern.hpp"

#include <cmath>
#if NIHU_MEX_DEBUGGING
#include <mex.h>
#endif


namespace NiHu
{

/** \brief store-wrapper of a statically stored line quadrature */
template <unsigned order>
struct line_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<line_domain> const quadrature;
};

/** \brief definition of the statically stored line quadrature member */
template <unsigned order>
gaussian_quadrature<line_domain> const line_quad_store<order>::quadrature(order);


/** \brief definition of Laurent coefficients of singularities
 * \tparam singuarity_type the singularity type
 * \details The class should implement a static function template eval templated on a guiggiani class
 * and taking a guiggiani object by reference as parameter.
 * The function template should compute the Laurent coefficients of the singularity in terms of
 * the Taylor series expansion of
 * - the distance vector
 * - the Jacobian vector
 * - the shape function vector
 *
 * and the unit normal at the singular point.
 */
template <class singularity_type>
class polar_laurent_coeffs;

/** \brief shorthand for constant minus two */
typedef std::integral_constant<int, -2> _m2;
/** \brief shorthand for constant minus one */
typedef std::integral_constant<int, -1> _m1;
/** \brief shorthand for constant zero */
typedef std::integral_constant<int, 0> _0;
/** \brief shorthand for constant one */
typedef std::integral_constant<int, 1> _1;
/** \brief shorthand for constant two */
typedef std::integral_constant<int, 2> _2;

/** \brief Implementation of Guiggiani's method
 * \tparam TrialField the trial field type
 * \tparam Kernel the kernel type
 * \tparam RadialOrder the quadrature order of radial integration
 * \tparam TangentialOrder the quadrature order of tangential integration
 */
template <class TrialField, class Kernel, unsigned RadialOrder, unsigned TangentialOrder = RadialOrder, class Enable = void>
class guiggiani;


/** brief specialisation of Guiggiani's method for surface elements */
template <class TrialField, class Kernel, unsigned RadialOrder, unsigned TangentialOrder>
class guiggiani<TrialField, Kernel, RadialOrder, TangentialOrder, typename std::enable_if<
	element_traits::is_surface_element<typename TrialField::elem_t>::value
>::type>
{
public:
	/** \brief quadrature orders stored as internal constants */
	enum {
		radial_order = RadialOrder,	/**< \brief quadrature order in radial direction */
		tangential_order = TangentialOrder	/**< \brief quadrature order in tangential direction */
	};

	/** \brief the trial field type */
	typedef TrialField trial_field_t;
	/** \brief the element type */
	typedef typename trial_field_t::elem_t elem_t;

	/** \brief the original reference domain type */
	typedef typename elem_t::domain_t domain_t;
	/** \brief the reference coordinate vector type */
	typedef typename domain_t::xi_t xi_t;
	/** \brief the geometrical scalar type */
	typedef typename elem_t::scalar_t scalar_t;
	/** \brief the Rong's transformation matrix type */
	typedef Eigen::Matrix<scalar_t, domain_t::dimension, domain_t::dimension> trans_t;

	/** \brief the physical coordinate vector type */
	typedef typename elem_t::x_t x_t;

	/** \brief shape function set type */
	typedef typename trial_field_t::nset_t trial_nset_t;
	/** \brief the shape function vector type */
	typedef typename trial_nset_t::shape_t trial_n_shape_t;

	/** \brief the kernel type */
	typedef Kernel kernel_t;

	/** \brief the kernel's test input type */
	typedef typename kernel_traits<kernel_t>::test_input_t test_input_t;
	/** \brief the kernel's trial input type */
	typedef typename kernel_traits<kernel_t>::trial_input_t trial_input_t;

	/** \brief the kernel's weighted trial input type */
	typedef typename merge<
		trial_input_t,
		typename build<normal_jacobian<typename trial_input_t::space_t> >::type
	>::type w_trial_input_t;

	/** \brief the singular core type */
	typedef typename singular_kernel_traits<kernel_t>::singular_core_t singular_core_t;
	/** \brief the Laurent coefficients computing class */
	typedef polar_laurent_coeffs<singular_core_t> laurent_t;

	/** \brief value type of the Laurent coefficients */
	typedef typename semi_block_product_result_type<
		typename singular_core_t::result_t, trial_n_shape_t
	>::type laurent_coeff_t;

	/** \brief value type of the integral result */
	typedef typename semi_block_product_result_type<
		typename kernel_t::result_t, trial_n_shape_t
	>::type total_result_t;

	/** \todo these calculations are only valid for the 3D case */
	enum {
		/** \brief the required Laurent expansion order */
		laurent_order = singular_kernel_traits<kernel_t>::singularity_type_t::value - 1,
	};

	/** \brief constructor
	 * \param [in] elem the element
	 * \param [in] kernel the kernel
	 */
	guiggiani(elem_t const &elem, kernel_base<kernel_t> const &kernel)
		: m_elem(elem), m_kernel(kernel)
	{
	}

private:
	/** \brief store xi0-related quantities for acceleration
	 * \param [in] xi0 the singular point in the reference domain
	 * \param [in] normal the normal in the singular point
	 */
	void compute_xi0(xi_t const &xi0, x_t const &normal)
	{
		// store the singular point and the unit normal vector
		m_xi0 = xi0;
		m_n0 = normal.normalized();

		// get the physical location and its gradient
		m_x0 = m_elem.get_x(m_xi0);
		typename elem_t::dx_return_type dx = m_elem.get_dx(m_xi0);

		scalar_t cosgamma = dx.col(shape_derivative_index::dXI).dot(dx.col(shape_derivative_index::dETA)) /
			dx.col(shape_derivative_index::dXI).norm() / dx.col(shape_derivative_index::dETA).norm();
		scalar_t u = dx.col(shape_derivative_index::dETA).norm() / dx.col(shape_derivative_index::dXI).norm();
		m_T <<
			1.0, cosgamma * u,
			0.0, std::sqrt(1.0-cosgamma*cosgamma) * u;
		// this scaling ensures that A is always equal to 1.0
		m_T *= dx.col(shape_derivative_index::dXI).norm();
		m_Tinv = m_T.inverse();

		m_eta0 = m_T * m_xi0;

		m_Jvec_series[0] = m_elem.get_normal(m_xi0) / m_T.determinant();
		m_N_series[0] = trial_nset_t::template eval_shape<0>(xi0);

		// geometrical parameters (planar helpers)
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			xi_t c1 = m_T * domain_t::get_corner(n);			// corner
			xi_t c2 = m_T * domain_t::get_corner((n + 1) % N);	// next corner
			xi_t l = xi_t::Zero();
			l(0) = 1.0;
			if ((c2-c1).norm() > 1e-3)
				l = (c2 - c1).normalized();					// side unit vector

			xi_t d1 = c1 - m_eta0;			// vector to corners
			xi_t d0 = d1 - l*d1.dot(l);		// perpendicular to side

			m_theta_lim[n] = std::atan2(d1(1), d1(0));	// corner angle
			m_theta0[n] = std::atan2(d0(1), d0(0));		// mid angle
			m_ref_distance[n] = d0.norm();				// distance to side
		}
	}

	/** \brief store theta-related quantities for acceleration
	 * \param [in] theta the angle in the reference domain
	 */
	void compute_theta(scalar_t theta)
	{
		// contains cos theta sin theta in xi system
		xi_t xi = m_Tinv.col(shape_derivative_index::dXI) * std::cos(theta)
			+ m_Tinv.col(shape_derivative_index::dETA) * std::sin(theta);

		typename elem_t::dx_return_type dx = m_elem.get_dx(m_xi0);

		m_rvec_series[0] = dx.col(shape_derivative_index::dXI) * xi(shape_derivative_index::dXI, 0)
			+ dx.col(shape_derivative_index::dETA) * xi(shape_derivative_index::dETA, 0);
		if (laurent_order > 1) // compile_time IF
		{
			typename elem_t::ddx_return_type ddx = m_elem.get_ddx(m_xi0);

			m_rvec_series[1] = ddx.col(shape_derivative_index::dXIXI) * xi(shape_derivative_index::dXI, 0)*xi(shape_derivative_index::dXI, 0) / 2.0 +
				ddx.col(shape_derivative_index::dXIETA) * xi(shape_derivative_index::dXI, 0)*xi(shape_derivative_index::dETA, 0) +
				ddx.col(shape_derivative_index::dETAETA) * xi(shape_derivative_index::dETA, 0)*xi(shape_derivative_index::dETA, 0) / 2.0;
				
			m_Jvec_series[1] = (
				xi(shape_derivative_index::dXI, 0) *
				(
				ddx.col(shape_derivative_index::dXIXI).cross(dx.col(shape_derivative_index::dETA)) +
				dx.col(shape_derivative_index::dXI).cross(ddx.col(shape_derivative_index::dXIETA))
				)
				+
				xi(shape_derivative_index::dETA, 0) *
				(
				ddx.col(shape_derivative_index::dETAXI).cross(dx.col(shape_derivative_index::dETA)) +
				dx.col(shape_derivative_index::dXI).cross(ddx.col(shape_derivative_index::dETAETA))
				)
				) / m_T.determinant();

			auto dN = trial_nset_t::template eval_shape<1>(m_xi0);
			m_N_series[1] = dN.col(shape_derivative_index::dXI)*xi(shape_derivative_index::dXI, 0)
				+ dN.col(shape_derivative_index::dETA)*xi(shape_derivative_index::dETA, 0);
		}
	}

public:
	/** \brief the entry point of integration
	 * \tparam result_t the result type where the integral is assembled
	 * \param [in] xi0 the singular point in intrinsic coordinates
	 * \param [in] normal the normal vector in the singular point
	 * \param [out] I the result matrix where the data is assembled
	 */
	template <class result_t>
	void integrate(result_t &&I, xi_t const &xi0, x_t const &normal)
	{
#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("Computing integral with Guiggiani.\nRadial_order: %d\nTangential order: %d\n", RadialOrder, TangentialOrder);
			printed = true;
		}
#endif

		// compute the xi0-related quantities
		compute_xi0(xi0, normal);

		// bind the kernel at the test input
		test_input_t test_input(m_elem, m_xi0);
		auto bound = m_kernel.bind(test_input);

		// iterate through triangles
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			// get angular integration limits
			scalar_t t1 = m_theta_lim[n];
			scalar_t t2 = m_theta_lim[(n + 1) % N];
			
			/** todo this is for highly distorted elements */
			if (std::abs(t2-t1) < 1e-3)
				continue;
			
			if (t2 < t1)
				t2 += 2.0 * M_PI;

			// theta integration
			for (auto const &q_theta : line_quad_store<tangential_order>::quadrature)
			{
				// compute theta base point and weight
				scalar_t xx = q_theta.get_xi()(0);
				scalar_t theta = ((1.0 - xx) * t1 + (1.0 + xx) * t2) / 2.0;
				scalar_t w_theta = q_theta.get_w() * (t2 - t1) / 2.0;

				// compute theta-related members
				compute_theta(theta);
				laurent_t::eval(*this);

				// reference domain's limit
				scalar_t rho_lim = m_ref_distance[n] / std::cos(theta - m_theta0[n]);

				laurent_coeff_t toadd = m_Fcoeffs[0] * std::log(rho_lim);
				if (laurent_order > 1) // compile_time IF
					toadd -= m_Fcoeffs[1] / rho_lim;
				toadd *= w_theta;

				for (int r = 0; r < I.rows(); ++r)	// loop needed for scalar casting
					for (int c = 0; c < I.cols(); ++c)
						I(r,c) += toadd(r,c);

				// radial part of surface integration
				for (auto const &q_rho : line_quad_store<radial_order>::quadrature)
				{
					// quadrature base point and weight
					scalar_t rho = (1.0 + q_rho.get_xi()(0)) * rho_lim / 2.0;
					scalar_t w_rho = q_rho.get_w() * rho_lim / 2.0;

					// compute location in original reference domain
					xi_t eta(rho*std::cos(theta), rho*std::sin(theta));
					xi_t xi = m_Tinv * eta + m_xi0;

					// evaluate G * N * J * rho
					w_trial_input_t trial_input(m_elem, xi);
					typename kernel_t::result_t GJrho =
						bound(trial_input) * (trial_input.get_jacobian() / m_T.determinant() * rho);
					typename trial_nset_t::shape_t N =
						trial_nset_t::template eval_shape<0>(xi);
					total_result_t F = semi_block_product(GJrho, N);

					// subtract the analytical singularity
					laurent_coeff_t singular_part = m_Fcoeffs[0];
					if (laurent_order > 1) // compile_time IF
						singular_part += m_Fcoeffs[1] / rho;
					singular_part /= rho;

					for (int r = 0; r < F.rows(); ++r)	// loop needed for scalar casting
						for (int c = 0; c < F.cols(); ++c)
							F(r,c) -= singular_part(r,c);

					// surface integral accumulation
					I += w_theta * w_rho * F;
				} // rho loop
			} // theta loop
		} // element sides
	} // end of function integrate

	/** \brief return Taylor coefficient of the distance measured from the collocation point */
	template <class T, T order>
	x_t const &get_rvec_series(std::integral_constant<T, order>) const
	{
		static_assert((order-1)>=0, "Required distance Taylor coefficient too low");
		static_assert((order-1) < laurent_order, "Required distance Taylor coefficient too high");
		return m_rvec_series[order-1];
	}

	/** \brief return Taylor coefficient of the Jacobian vector around the collocation point */
	template <class T, T order>
	x_t const &get_Jvec_series(std::integral_constant<T, order>) const
	{
		static_assert(order < laurent_order, "Required Jacobian Taylor coefficient too high");
		static_assert(order >= 0, "Required Jacobian Taylor coefficient too low");
		return m_Jvec_series[order];
	}

	/** \brief return Taylor coefficient of the shape set around the collocation point */
	template <class T, T order>
	trial_n_shape_t const &get_shape_series(std::integral_constant<T, order>) const
	{
		static_assert(order < laurent_order, "Required shape set Taylor coefficient too high");
		static_assert(order >= 0, "Required Jacobian Taylor coefficient too low");
		return m_N_series[order];
	}

	/** \brief return a Laurent coefficient */
	template <class T, T order>
	laurent_coeff_t &get_laurent_coeff(std::integral_constant<T, order>)
	{
		static_assert(-(order+1) < laurent_order, "Required Laurent coefficient too high");
		static_assert(-(order+1) >= 0, "Required Laurent coefficient too low");
		return m_Fcoeffs[-(order+1)];
	}

	/** \brief set a Laurent coefficient */
	template <class T, T order>
	void set_laurent_coeff(std::integral_constant<T, order>, laurent_coeff_t const &v)
	{
		static_assert(-(order+1) < laurent_order, "Required Laurent coefficient too high");
		static_assert(-(order+1) >= 0, "Required Laurent coefficient too low");
		m_Fcoeffs[-(order+1)] = v;
	}

	/** \brief return the unit normal at the collocation point */
	x_t const &get_n0(void) const
	{
		return m_n0;
	}
	
	kernel_t const &get_kernel(void) const
	{
		return m_kernel.derived();
	}
	
private:
	elem_t const &m_elem;		/**< \brief the element reference */
	kernel_base<kernel_t> const &m_kernel;	/**< \brief the kernel reference */

	xi_t m_xi0;		/**< \brief the source local coordinate */
	x_t m_x0;		/**< \brief the source point */
	x_t m_n0;		/**< \brief the unit normal vector at the source point */

	trans_t m_T;	/**< \brief transformation to distorted reference domain */
	trans_t m_Tinv;	/**< \brief inverse of the transformation */
	xi_t m_eta0;	/**< \brief the source local coordinate in the distorted reference domain */

	/** \brief angles to the corners of the distorted reference domain */
	scalar_t m_theta_lim[domain_t::num_corners];
	/** \brief angles to the sides of the distorted reference domain */
	scalar_t m_theta0[domain_t::num_corners];
	/** \brief distances to the distorted reference domain */
	scalar_t m_ref_distance[domain_t::num_corners];

	x_t m_rvec_series[laurent_order];	        /**< \brief series expansion of the location vector */
	x_t m_Jvec_series[laurent_order];	        /**< \brief series expansion of the Jacobian vector */
	trial_n_shape_t m_N_series[laurent_order];	/**< \brief series expansion of the the shape function vector */
	laurent_coeff_t m_Fcoeffs[laurent_order];	/**< \brief the Laurent coefficient */
};



template <class TrialField, class Kernel, unsigned RadialOrder, unsigned TangentialOrder>
class guiggiani<TrialField, Kernel, RadialOrder, TangentialOrder, typename std::enable_if<
	!element_traits::is_surface_element<typename TrialField::elem_t>::value
>::type>
{
public:
	/** \brief quadrature orders stored as internal constants */
	enum {
		radial_order = RadialOrder,	/**< \brief quadrature order in radial direction */
		tangential_order = TangentialOrder	/**< \brief quadrature order in tangential direction */
	};

	/** \brief the trial field type */
	typedef TrialField trial_field_t;
	/** \brief the element type */
	typedef typename trial_field_t::elem_t elem_t;

	/** \brief the original reference domain type */
	typedef typename elem_t::domain_t domain_t;
	/** \brief the reference coordinate vector type */
	typedef typename domain_t::xi_t xi_t;
	/** \brief the geometrical scalar type */
	typedef typename elem_t::scalar_t scalar_t;
	/** \brief the Rong's transformation matrix type */
	typedef Eigen::Matrix<scalar_t, domain_t::dimension, domain_t::dimension> trans_t;

	/** \brief the physical coordinate vector type */
	typedef typename elem_t::x_t x_t;

	/** \brief shape function set type */
	typedef typename trial_field_t::nset_t trial_nset_t;
	/** \brief the shape function vector type */
	typedef typename trial_nset_t::shape_t trial_n_shape_t;

	/** \brief the kernel type */
	typedef Kernel kernel_t;

	/** \brief the kernel's test input type */
	typedef typename kernel_traits<kernel_t>::test_input_t test_input_t;
	/** \brief the kernel's trial input type */
	typedef typename kernel_traits<kernel_t>::trial_input_t trial_input_t;

	/** \brief the kernel's weighted trial input type */
	typedef typename merge<
		trial_input_t,
		typename build<volume_jacobian<typename trial_input_t::space_t> >::type
	>::type w_trial_input_t;

	/** \brief the singular core type */
	typedef typename singular_kernel_traits<kernel_t>::singular_core_t singular_core_t;
	/** \brief the Laurent coefficients computing class */
	typedef polar_laurent_coeffs<singular_core_t> laurent_t;

	/** \brief value type of the Laurent coefficients */
	typedef typename semi_block_product_result_type<
		typename singular_core_t::result_t, trial_n_shape_t
	>::type laurent_coeff_t;

	/** \brief value type of the integral result */
	typedef typename semi_block_product_result_type<
		typename kernel_t::result_t, trial_n_shape_t
	>::type total_result_t;

	/** \todo these calculations are only valid for the 2D volume case */
	enum {
		/** \brief the required Laurent expansion order */
		laurent_order = singular_kernel_traits<kernel_t>::singularity_type_t::value - 1,
	};

	/** \brief constructor
	 * \param [in] elem the element
	 * \param [in] kernel the kernel
	 */
	guiggiani(elem_t const &elem, kernel_base<kernel_t> const &kernel)
		: m_elem(elem), m_kernel(kernel)
	{
	}

private:
	/** \brief store xi0-related quantities for acceleration
	 * \param [in] xi0 the singular point in the reference domain
	 */
	void compute_xi0(xi_t const &xi0)
	{
		m_xi0 = xi0;
		m_x0 = m_elem.get_x(m_xi0);
		typename elem_t::dx_return_type dx = m_elem.get_dx(m_xi0);

		scalar_t cosgamma = dx.col(shape_derivative_index::dXI).dot(dx.col(shape_derivative_index::dETA)) /
			dx.col(shape_derivative_index::dXI).norm() / dx.col(shape_derivative_index::dETA).norm();
		scalar_t u = dx.col(shape_derivative_index::dETA).norm() / dx.col(shape_derivative_index::dXI).norm();
		m_T <<
			1.0, cosgamma * u,
			0.0, std::sqrt(1.0-cosgamma*cosgamma) * u;
		// this scaling ensures that A is always equal to 1.0
		m_T *= dx.col(shape_derivative_index::dXI).norm();
		m_Tinv = m_T.inverse();

		m_eta0 = m_T * m_xi0;

		m_J_series[0] = m_elem.get_dx(m_xi0).determinant() / m_T.determinant();
		m_N_series[0] = trial_nset_t::template eval_shape<0>(xi0);

		// geometrical parameters (planar helpers)
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			xi_t c1 = m_T * domain_t::get_corner(n);			// corner
			xi_t c2 = m_T * domain_t::get_corner((n + 1) % N);	// next corner
			xi_t l = (c2 - c1).normalized();					// side vector

			xi_t d1 = c1 - m_eta0;			// vector to corners
			xi_t d0 = d1 - l*d1.dot(l);		// perpendicular to side

			m_theta_lim[n] = std::atan2(d1(1), d1(0));	// corner angle
			m_theta0[n] = std::atan2(d0(1), d0(0));		// mid angle
			m_ref_distance[n] = d0.norm();				// distance to side
		}
	}

	template <class Input1, class Input2>
	typename Input1::Scalar detvectors(Eigen::DenseBase<Input1> const &v1, Eigen::DenseBase<Input2> const &v2)
	{
		return v1(0)*v2(1) - v1(1)*v2(0);
	}

	/** \brief store theta-related quantities for acceleration
	 * \param [in] theta the angle in the reference domain
	 */
	void compute_theta(scalar_t theta)
	{
		// contains cos theta sin theta in xi system
		xi_t xi = m_Tinv.col(shape_derivative_index::dXI) * std::cos(theta)
			+ m_Tinv.col(shape_derivative_index::dETA) * std::sin(theta);

		typename elem_t::dx_return_type dx = m_elem.get_dx(m_xi0);

		m_rvec_series[0] = dx.col(shape_derivative_index::dXI) * xi(shape_derivative_index::dXI, 0)
			+ dx.col(shape_derivative_index::dETA) * xi(shape_derivative_index::dETA, 0);
		if (laurent_order > 1) // compile_time IF
		{
			typename elem_t::ddx_return_type ddx = m_elem.get_ddx(m_xi0);

			m_rvec_series[1] =
				ddx.col(shape_derivative_index::dXIXI)*
				xi(shape_derivative_index::dXI, 0)*xi(shape_derivative_index::dXI, 0) / 2.0
				+
				ddx.col(shape_derivative_index::dXIETA)*
				xi(shape_derivative_index::dXI, 0)*xi(shape_derivative_index::dETA, 0)
				+
				ddx.col(shape_derivative_index::dETAETA)*
				xi(shape_derivative_index::dETA, 0)*xi(shape_derivative_index::dETA, 0) / 2.0;

			m_J_series[1] = (
				xi(0, 0) *
				(
				detvectors(ddx.col(shape_derivative_index::dXIXI), dx.col(shape_derivative_index::dETA)) +
				detvectors(dx.col(shape_derivative_index::dXI), ddx.col(shape_derivative_index::dXIETA))
				)
				+
				xi(1, 0) *
				(
				detvectors(ddx.col(shape_derivative_index::dETAXI), dx.col(shape_derivative_index::dETA))+
				detvectors(dx.col(shape_derivative_index::dXI), ddx.col(shape_derivative_index::dETAETA))
				)
				) / m_T.determinant();

			auto dN = trial_nset_t::template eval_shape<1>(m_xi0);
			m_N_series[1] = dN.col(shape_derivative_index::dXI)*xi(shape_derivative_index::dXI, 0)
				+ dN.col(shape_derivative_index::dETA)*xi(shape_derivative_index::dETA, 0);
		}
	}

public:
	/** \brief the entry point of integration
	 * \tparam result_t the result type where the integral is assembled
	 * \param [in] xi0 the singular point in intrinsic coordinates
	 * \param [in] normal the normal vector in the singular point
	 * \param [out] I the result matrix where the data is assembled
	 */
	template <class result_t>
	void integrate(result_t &&I, xi_t const &xi0)
	{
#if NIHU_MEX_DEBUGGING
		static bool printed = false;
		if (!printed)
		{
			mexPrintf("Computing integral with Guiggiani.\nRadial_order: %d\nTangential order: %d\n", RadialOrder, TangentialOrder);
			printed = true;
		}
#endif
		
		// compute the xi0-related quantities
		compute_xi0(xi0);

		// bind the kernel at the test input
		test_input_t test_input(m_elem, m_xi0);
		auto bound = m_kernel.bind(test_input);

		// iterate through triangles
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			// get angular integration limits
			scalar_t t1 = m_theta_lim[n];
			scalar_t t2 = m_theta_lim[(n + 1) % N];
			if (t2 < t1)
				t2 += 2.0 * M_PI;

			// theta integration
			for (auto const &q_theta : line_quad_store<tangential_order>::quadrature)
			{
				// compute theta base point and weight
				scalar_t xx = q_theta.get_xi()(0);
				scalar_t theta = ((1.0 - xx) * t1 + (1.0 + xx) * t2) / 2.0;
				scalar_t w_theta = q_theta.get_w() * (t2 - t1) / 2.0;

				// compute theta-related members
				compute_theta(theta);
				laurent_t::eval(*this);

				// reference domain's limit
				scalar_t rho_lim = m_ref_distance[n] / std::cos(theta - m_theta0[n]);

				laurent_coeff_t toadd = m_Fcoeffs[0] * std::log(rho_lim);
				if (laurent_order > 1) // compile_time IF
					toadd -= m_Fcoeffs[1] / rho_lim;
				toadd *= w_theta;

				for (int r = 0; r < I.rows(); ++r)	// loop needed for scalar casting
					for (int c = 0; c < I.cols(); ++c)
						I(r,c) += toadd(r,c);

				// radial part of surface integration
				for (auto const &q_rho : line_quad_store<radial_order>::quadrature)
				{
					// quadrature base point and weight
					scalar_t rho = (1.0 + q_rho.get_xi()(0)) * rho_lim / 2.0;
					scalar_t w_rho = q_rho.get_w() * rho_lim / 2.0;

					// compute location in original reference domain
					xi_t eta(rho*std::cos(theta), rho*std::sin(theta));
					xi_t xi = m_Tinv * eta + m_xi0;

					// evaluate G * N * J * rho
					w_trial_input_t trial_input(m_elem, xi);
					typename kernel_t::result_t GJrho =
						bound(trial_input) * (trial_input.get_jacobian() / m_T.determinant() * rho);
					typename trial_nset_t::shape_t N =
						trial_nset_t::template eval_shape<0>(xi);
					total_result_t F = semi_block_product(GJrho, N);

					// subtract the analytical singularity
					laurent_coeff_t singular_part = m_Fcoeffs[0];
					if (laurent_order > 1) // compile_time IF
						singular_part += m_Fcoeffs[1] / rho;
					singular_part /= rho;

					for (int r = 0; r < F.rows(); ++r)	// loop needed for scalar casting
						for (int c = 0; c < F.cols(); ++c)
							F(r,c) -= singular_part(r,c);

					// surface integral accumulation
					I += w_theta * w_rho * F;
				} // rho loop
			} // theta loop
		} // element sides
	}

	/** \brief return Taylor coefficient of the distance measured from the collocation point */
	template <class T, T order>
	x_t const &get_rvec_series(std::integral_constant<T, order>) const
	{
		static_assert((order-1)>=0, "Required distance Taylor coefficient too low");
		static_assert((order-1) < laurent_order, "Required distance Taylor coefficient too high");
		return m_rvec_series[order-1];
	}

	/** \brief return Taylor coefficient of the Jacobian vector around the collocation point */
	template <class T, T order>
	double get_J_series(std::integral_constant<T, order>) const
	{
		static_assert(order < laurent_order, "Required Jacobian Taylor coefficient too high");
		static_assert(order >= 0, "Required Jacobian Taylor coefficient too low");
		return m_J_series[order];
	}

	/** \brief return Taylor coefficient of the shape set around the collocation point */
	template <class T, T order>
	trial_n_shape_t const &get_shape_series(std::integral_constant<T, order>) const
	{
		static_assert(order < laurent_order, "Required shape set Taylor coefficient too high");
		static_assert(order >= 0, "Required Jacobian Taylor coefficient too low");
		return m_N_series[order];
	}

	/** \brief return a Laurent coefficient */
	template <class T, T order>
	laurent_coeff_t &get_laurent_coeff(std::integral_constant<T, order>)
	{
		static_assert(-(order+1) < laurent_order, "Required Laurent coefficient too high");
		static_assert(-(order+1) >= 0, "Required Laurent coefficient too low");
		return m_Fcoeffs[-(order+1)];
	}

	/** \brief set a Laurent coefficient */
	template <class T, T order>
	void set_laurent_coeff(std::integral_constant<T, order>, laurent_coeff_t const &v)
	{
		static_assert(-(order+1) < laurent_order, "Required Laurent coefficient too high");
		static_assert(-(order+1) >= 0, "Required Laurent coefficient too low");
		m_Fcoeffs[-(order+1)] = v;
	}
	
	
	kernel_t const &get_kernel(void) const
	{
		return m_kernel.derived();
	}
	

private:
	elem_t const &m_elem;		/**< \brief the element reference */
	kernel_base<kernel_t> const &m_kernel;	/**< \brief the kernel reference */

	xi_t m_xi0;		/**< \brief the source local coordinate */
	x_t m_x0;		/**< \brief the source point */

	trans_t m_T;	/**< \brief transformation to distorted reference domain */
	trans_t m_Tinv;	/**< \brief inverse of the transformation */
	xi_t m_eta0;	/**< \brief the source local coordinate in the distorted reference domain */

	/** \brief angles to the corners of the distorted reference domain */
	scalar_t m_theta_lim[domain_t::num_corners];
	/** \brief angles to the sides of the distorted reference domain */
	scalar_t m_theta0[domain_t::num_corners];
	/** \brief distances to the distorted reference domain */
	scalar_t m_ref_distance[domain_t::num_corners];

	x_t m_rvec_series[laurent_order];	        /**< \brief series expansion of the location vector */
	scalar_t m_J_series[laurent_order];	        	/**< \brief series expansion of the Jacobian vector */
	trial_n_shape_t m_N_series[laurent_order];	/**< \brief series expansion of the the shape function vector */
	laurent_coeff_t m_Fcoeffs[laurent_order];	/**< \brief the Laurent coefficient */
};

}



#endif // GUIGGIANI_1992_HPP_INCLUDED

