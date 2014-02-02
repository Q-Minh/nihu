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

/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for hypersingular collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include <cmath>
#include "../util/store_pattern.hpp"
#include "../core/field.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "location_normal.hpp"

/** \brief store-wrapper of a statically stored quadrature */
template <unsigned order>
struct line_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<line_domain> const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <unsigned order>
gaussian_quadrature<line_domain> const line_quad_store<order>::quadrature(order);


/** \brief definition of Laurent coefficients of singularities
 * \tparam singuarity_type the singularity type
 * \details The class should implement a static function template eval templated on a guiggiani class
 * and taking a guiggiani object by reference as parameter.
 * The function  template should compute the Laurent coefficients of the singularity in terms of
 * the series expansion of
 * - the distance vector
 * - the Jacobian vector
 * - the shape function vector
 *
 * and the unit normal at the singular point.
 */
template <class singularity_type>
class polar_laurent_coeffs;


/** \brief Implementation of Guiggiani's method
 * \tparam TrialField the trial field type
 * \tparam Kernel the kernel type
 * \tparam RadialOrder the quadrature order of radial integration
 * \tparam TangentialOrder the quadrature order of tangential integration
 */
template <class TrialField, class Kernel, unsigned RadialOrder, unsigned TangentialOrder = RadialOrder>
class guiggiani
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

	/** \brief the singular kernel ancestor type */
	typedef typename singular_kernel_traits<kernel_t>::singular_kernel_ancestor_t singular_kernel_ancestor_t;
	/** \brief the Laurent coefficients computing class */
	typedef polar_laurent_coeffs<singular_kernel_ancestor_t> laurent_t;

	template <class singularity_type>
	friend class polar_laurent_coeffs;

	/** \brief constructor
	 * \param [in] elem the element
	 * \param [in] kernel the kernel
	 */
	guiggiani(elem_t const &elem, kernel_t const &kernel)
		: m_elem(elem), m_kernel(kernel)
	{
	}

private:
	/** \brief store xi0-related quantities for acceleration
	 * \param [in] xi0 the singular point in the reference domain
	 */
	void compute_xi0(xi_t const &xi0, x_t const &normal)
	{
		m_xi0 = xi0;
		m_n0 = normal.normalized();

		m_x0 = m_elem.get_x(m_xi0);
		typename elem_t::dx_t dx = m_elem.get_dx(m_xi0);

		scalar_t cosgamma = dx.col(shape_index::dXI).dot(dx.col(shape_index::dETA)) / dx.col(shape_index::dXI).norm() / dx.col(shape_index::dETA).norm();
		scalar_t u = dx.col(shape_index::dETA).norm() / dx.col(shape_index::dXI).norm();
		m_T <<
			1.0, cosgamma * u,
			0.0, std::sqrt(1.0-cosgamma*cosgamma) * u;
		m_Tinv = m_T.inverse();

		m_eta0 = m_T * m_xi0;

		m_Jvec_series[0] = m_elem.get_normal(m_xi0) / m_T.determinant();
		m_N_series[0] = trial_nset_t::eval_shape(xi0);

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

	/** \brief store theta-related quantities for acceleration
	 * \param [in] theta the angle in the reference domain
	 */
	void compute_theta(scalar_t theta)
	{
		// contains cos theta sin theta in xi system
		xi_t xi = m_Tinv.col(shape_index::dXI) * std::cos(theta) + m_Tinv.col(shape_index::dETA) * std::sin(theta);

		typename elem_t::dx_t dx = m_elem.get_dx(m_xi0);
		typename elem_t::ddx_t ddx = m_elem.get_ddx(m_xi0);

		m_rvec_series[0] = dx.col(shape_index::dXI) * xi(shape_index::dXI) + dx.col(shape_index::dETA) * xi(shape_index::dETA);
		m_A = m_rvec_series[0].norm();
		m_rvec_series[1] = ddx.col(shape_index::dXIXI) * xi(shape_index::dXI)*xi(shape_index::dXI) / 2.0 +
			ddx.col(shape_index::dXIETA) * xi(shape_index::dXI)*xi(1) +
			ddx.col(shape_index::dETAETA) * xi(shape_index::dETA)*xi(shape_index::dETA) / 2.0;

		m_Jvec_series[1] = (
			xi(shape_index::dXI) * (ddx.col(shape_index::dXIXI).cross(dx.col(shape_index::dETA)) + dx.col(shape_index::dXI).cross(ddx.col(shape_index::dXIETA))) +
			xi(shape_index::dETA) * (ddx.col(shape_index::dETAXI).cross(dx.col(shape_index::dETA)) + dx.col(shape_index::dXI).cross(ddx.col(shape_index::dETAETA)))
			) / m_T.determinant();

		auto dN = trial_nset_t::eval_dshape(m_xi0);
		m_N_series[1] = dN.col(shape_index::dXI)*xi(shape_index::dXI) + dN.col(shape_index::dETA)*xi(shape_index::dETA);
	}

public:
	/** \brief the entry point of integration
	 * \tparam result_t the result type where the integral is assembled
	 * \param [in] xi0 the singular point in intrinsic coordinates
	 * \param [out] result the result matrix where the data is assembled
	 */
	template <class result_t>
	void integrate(result_t &I, xi_t const &xi0, x_t const &normal)
	{
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
			if (t2 < t1)
				t2 += 2.0 * M_PI;

			// theta integration
			auto const &quad_theta = line_quad_store<tangential_order>::quadrature;
			for (auto it_theta = quad_theta.begin(); it_theta != quad_theta.end(); ++it_theta)
			{
				// compute theta base point and weight
				scalar_t xx = it_theta->get_xi()(0);
				scalar_t theta = ((1.0 - xx) * t1 + (1.0 + xx) * t2) / 2.0;
				scalar_t w_theta = it_theta->get_w() * (t2 - t1) / 2.0;

				// compute theta-related members
				compute_theta(theta);
				laurent_t::eval(*this);

				// reference domain's limit
				scalar_t rho_lim = m_ref_distance[n] / std::cos(theta - m_theta0[n]);

				auto toadd = w_theta * (m_Fcoeffs[0] * std::log(rho_lim) - m_Fcoeffs[1] * (1.0 / rho_lim));
				for (int j = 0; j < I.cols(); ++j)	// loop needed for scalar casting
					I(j) += toadd(j);

				// radial part of surface integration
				auto const &quad_rho = line_quad_store<radial_order>::quadrature;
				for (auto it_rho = quad_rho.begin(); it_rho != quad_rho.end(); ++it_rho)
				{
					// quadrature base point and weight
					scalar_t rho = (1.0 + it_rho->get_xi()(0)) * rho_lim / 2.0;
					scalar_t w_rho = it_rho->get_w() * rho_lim / 2.0;

					// compute location in original reference domain
					xi_t eta(rho*std::cos(theta), rho*std::sin(theta));
					xi_t xi = m_Tinv * eta + m_xi0;

					// evaluate G * N * J * rho
					w_trial_input_t trial_input(m_elem, xi);
					auto F = (
						bound(trial_input) *
						trial_input.get_jacobian() / m_T.determinant() *
						rho *
						trial_nset_t::eval_shape(xi)
						).eval();

					// subtract the analytical singularity
					auto singular_part = (m_Fcoeffs[1] / rho + m_Fcoeffs[0]) / rho;
					for (int j = 0; j < F.cols(); ++j)	// loop needed for scalar casting
						F(j) -= singular_part(j);

					// surface integral accumulation
					I += w_theta * w_rho * F;
				} // rho loop
			} // theta loop
		} // element sides
	}

protected:
	elem_t const &m_elem;	/**< \brief the element reference */
	kernel_t const &m_kernel;	/**< \brief the kernel reference */

	xi_t m_xi0;				/**< \brief the source local coordinate */
	x_t m_x0;				/**< \brief the source point */
	x_t m_n0;				/**< \brief the unit normal vector at the source point */

	trans_t m_T;			/**< \brief transformation to distorted reference domain */
	trans_t m_Tinv;			/**< \brief inverse of the transformation */
	xi_t m_eta0;			/**< \brief the source local coordinate in the distorted reference domain */

	/** \brief angles to the corners of the distorted reference domain */
	scalar_t m_theta_lim[domain_t::num_corners];
	/** \brief angles to the sides of the distorted reference domain */
	scalar_t m_theta0[domain_t::num_corners];
	/** \brief distances to the distorted reference domain */
	scalar_t m_ref_distance[domain_t::num_corners];

	scalar_t m_A;			        /**< \brief the magnitude of the location derivative */

	x_t m_rvec_series[2];	        /**< \brief series expansion of the location vector */
	x_t m_Jvec_series[2];	        /**< \brief series expansion of the Jacobian vector */
	trial_n_shape_t m_N_series[2];	/**< \brief series expansion of the the shape function vector */
	trial_n_shape_t m_Fcoeffs[2];	/**< \brief the 1st and 2nd order Laurent coefficient */
};

#include "../library/laplace_kernel.hpp"

/** \brief specialisation of class ::polar_laurent_coeffs for the ::laplace_3d_HSP_kernel */
template <>
class polar_laurent_coeffs<laplace_3d_HSP_kernel>
{
public:
	template <class guiggiani>
	static void eval(guiggiani &obj)
	{
        auto A2 = obj.m_A * obj.m_A, A3 = A2 * obj.m_A;
		auto g1vec = obj.m_rvec_series[0] / A2 * (obj.m_rvec_series[1].dot(obj.m_Jvec_series[0]) + obj.m_rvec_series[0].dot(obj.m_Jvec_series[1]));

		auto b0vec = -obj.m_Jvec_series[0];
		auto b1vec = 3.0 * g1vec - obj.m_Jvec_series[1];

		auto a0 = b0vec.dot(obj.m_n0) * obj.m_N_series[0];
		auto a1 = b1vec.dot(obj.m_n0) * obj.m_N_series[0] + b0vec.dot(obj.m_n0) * obj.m_N_series[1];

		auto Sm2 = -3.0 * obj.m_rvec_series[0].dot(obj.m_rvec_series[1]) / A2 / A3;
		auto Sm3 = 1.0 / A3;

		obj.m_Fcoeffs[0] = -(Sm2 * a0 + Sm3 * a1) / (4.0 * M_PI);
		obj.m_Fcoeffs[1] = -(Sm3 * a0) / (4.0 * M_PI);
	}
};

#endif // GUIGGIANI_1992_HPP_INCLUDED
