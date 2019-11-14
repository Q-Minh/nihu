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

/**
 * \file laplace_singular_integrals.hpp
 * \brief Analytical expressions for the singular integrals of Laplace kernels
 * \ingroup lib_laplace
 *
 * \details
 * Analytical expression for the Laplace kernels over plane elements.
 */
#ifndef LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED
#define LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../core/match_types.hpp"
#include "../core/singular_integral_shortcut.hpp"
#include "guiggiani_1992.hpp"
#include "lib_element.hpp"
#include "laplace_kernel.hpp"
#include "normal_derivative_singular_integrals.hpp"
#include "plane_element_helper.hpp"
#include "quadrature_store_helper.hpp"

#include <boost/math/constants/constants.hpp>

#include <iostream>

namespace NiHu
{

/** \brief Collocational integral of the 2D SLP kernel over a curved line with general shape set
 * The log type singularity is subtracted in the reference domain.
 * The singularity is integrated analytically.
 * The regular part is integrated using Gaussian quadratures
 */
template <class TestField, class TrialField, size_t order>
class laplace_2d_SLP_collocation_general
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;

	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;

	typedef Eigen::Matrix<double, nTest, nTrial> result_t;

	typedef typename trial_field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;

	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;

	typedef regular_quad_store<domain_t, order> quadrature_t;

public:
	static result_t eval(elem_t const &elem)
	{
		using namespace boost::math::double_constants;

		result_t result = result_t::Zero();

		xi_t const &a = domain_t::get_corner(0);
		xi_t const &b = domain_t::get_corner(1);

		// traverse collocation points
		for (size_t i = 0; i < nTest; ++i)
		{
			xi_t const &xi0 = test_shape_t::corner_at(i);

			// trial shape function at the singular point
			auto N0 = trial_shape_t::template eval_shape<0>(xi0);
			auto N1 = trial_shape_t::template eval_shape<1>(xi0);
			auto N2 = trial_shape_t::template eval_shape<2>(xi0) / 2.;

			// singular point and jacobians
			x_t x = elem.get_x(xi0);
			double jac0 = elem.get_dx(xi0).norm();
			double jac1 = elem.get_dx(xi0).dot(elem.get_ddx(xi0)) / jac0;

			auto C0 = N0 * jac0;
			auto C1 = (N1 * jac0 + N0 * jac1);
			auto C2 = (N2 * jac0 + N1 * jac1);

			// traverse quadrature points
			for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
			{
				// transform quadrature to [xi0 b];
				xi_t xi = it->get_xi() * ((b(0) - xi0(0)) / 2.) + (xi0 + b) / 2.;
				double w = it->get_w() * ((b(0) - xi0(0)) / 2.);

				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				double jac = elem.get_dx(xi).norm();

				// evaluate Green's function
				double G = laplace_2d_SLP_kernel()(x, y);
				// multiply Shape function
				auto N = trial_shape_t::template eval_shape<0>(xi);
				// evaluate integrand
				auto F = G * N * jac;
				// evaluate integrand's singular part
				double rho = xi(0) - xi0(0);
				auto F0 = -std::log(std::abs(rho) * jac0) / two_pi * (C0 + rho * (C1 + rho * C2));

				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}

			// traverse quadrature points
			for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
			{
				// transform quadrature to  [a xi0]
				xi_t xi = it->get_xi() * ((xi0(0) - a(0)) / 2.) + (xi0 + a) / 2.;
				double w = it->get_w() * ((xi0(0) - a(0)) / 2.);

				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				double jac = elem.get_dx(xi).norm();

				// evaluate Green's function
				double G = laplace_2d_SLP_kernel()(x, y);
				// Shape function
				auto N = trial_shape_t::template eval_shape<0>(xi);
				// evaluate integrand
				auto F = G * N * jac;
				// evaluate singular part
				double rho = xi(0) - xi0(0);
				auto F0 = -std::log(std::abs(rho) * jac0) / two_pi * (C0 + rho * (C1 + rho * C2));

				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}

			// add analytic integral of singular part
			double rho1 = std::abs(a(0) - xi0(0));
			double rho2 = std::abs(b(0) - xi0(0));
			double d1 = rho1 * jac0;
			double d2 = rho2 * jac0;
			result.row(i) += 1. / two_pi * (
				C0 * (rho2 * (1. - std::log(d2)) + rho1 * (1. - std::log(d1)))
				+
				C1 / 4. * (rho2 * rho2 * (1. - 2. * std::log(d2)) - rho1 * rho1 * (1. - 2. * std::log(d1)))
				+
				C2 / 9. * (rho2 * rho2 * rho2 * (1. - 3. * std::log(d2)) + rho1 * rho1 * rho1 * (1. - 3. * std::log(d1)))
				);
		} // loop over collocation points

		return result;
	}
};


/**
 \brief Collocational integral of the 2D SLP kernel over a straight line with second order shape sets
 This is the simplification of function laplace_2d_SLP_collocation_general, skipping the Gaussian quadrature
 */
template <class TestField, class TrialField>
class laplace_2d_SLP_collocation_straight_line_second_order
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	typedef typename test_field_t::nset_t test_shape_set_t;
	typedef typename trial_field_t::nset_t trial_shape_set_t;

	typedef line_1_elem elem_t;

	typedef typename test_shape_set_t::xi_t xi_t;

	static size_t const rows = test_shape_set_t::num_nodes;
	static size_t const cols = trial_shape_set_t::num_nodes;
	static size_t const order = trial_shape_set_t::polynomial_order;

	typedef Eigen::Matrix<double, rows, cols> result_t;

public:
	static result_t eval(elem_t const &elem)
	{
		using namespace boost::math::double_constants;

		// get Jacobian
		double jac = elem.get_normal().norm();

		result_t result;

		// traverse collocational points
		for (size_t i = 0; i < rows; ++i)
		{
			xi_t const &xi0 = test_shape_set_t::corner_at(i);
			double rho1 = (1. + xi0(0));
			double rho2 = (1. - xi0(0));
			double d1 = jac * rho1;
			double d2 = jac * rho2;

			auto N0 = trial_shape_set_t::template eval_shape<0>(xi0);
			result.row(i) = N0 * (d2 * (1. - std::log(d2)) + d1 * (1. - std::log(d1)));

			if (trial_shape_set_t::polynomial_order >= 1)
			{
				auto N1 = trial_shape_set_t::template eval_shape<1>(xi0) / jac;
				result.row(i) += N1 * (d2 * d2 / 4. * (1. - 2. * std::log(d2)) - d1 * d1 / 4. * (1. - 2. * std::log(d1)));
			}

			if (trial_shape_set_t::polynomial_order >= 2)
			{
				auto N2 = trial_shape_set_t::template eval_shape<2>(xi0) / jac / jac / 2.;
				result.row(i) += N2 * (d2 * d2 * d2 / 9. * (1. - 3. * std::log(d2)) + d1 * d1 * d1 / 9. * (1. - 3. * std::log(d1)));
			}
		}

		return result / two_pi;
	}
};


/** \brief Collocational integral of the 2D HSP kernel over a general curved line with general shape sets */
template <class TestField, class TrialField, size_t order>
class laplace_2d_HSP_collocation_general
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;

	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;

	typedef Eigen::Matrix<double, nTest, nTrial> result_t;

	typedef typename trial_field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;

	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;

	typedef regular_quad_store<domain_t, order> quadrature_t;

public:

	static result_t eval(elem_t const &elem)
	{
		using namespace boost::math::double_constants;

		result_t result = result_t::Zero();

		xi_t const &a = domain_t::get_corner(0);
		xi_t const &b = domain_t::get_corner(1);

		// traverse collocation points
		for (size_t i = 0; i < nTest; ++i)
		{
			xi_t xi0 = test_shape_t::corner_at(i);

			// trial shape function and its derivative at singular point
			auto N0 = trial_shape_t::template eval_shape<0>(xi0);
			auto N1 = trial_shape_t::template eval_shape<1>(xi0);

			// singular point and normal at singular point
			x_t x = elem.get_x(xi0);
			x_t Jxvec = elem.get_normal(xi0);
			x_t nx = Jxvec.normalized();
			double jac0 = Jxvec.norm();
			double twopiJ0 = two_pi * jac0;

			// traverse quadrature points
			for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
			{
				// transform quadrature to [xi0 b];
				xi_t xi = it->get_xi() * ((b(0) - xi0(0)) / 2.) + (xi0 + b) / 2.;
				double w = it->get_w() * ((b(0) - xi0(0)) / 2.);
				double rho = xi(0) - xi0(0);

				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();

				// evaluate Green's function
				double G = laplace_2d_HSP_kernel()(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = G * N * jac;
				// evaluate singular part
				auto F0 = (N0 / rho + N1) / rho / twopiJ0;

				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}

			// traverse quadrature points
			for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
			{
				// transform quadrature to  [a xi0]
				xi_t xi = it->get_xi() * ((xi0(0) - a(0)) / 2.) + (xi0 + a) / 2.;
				double w = it->get_w() * ((xi0(0) - a(0)) / 2.);
				double rho = xi(0) - xi0(0);

				// get trial location, jacobian and normal
				x_t y = elem.get_x(xi);
				x_t Jyvec = elem.get_normal(xi);
				x_t ny = Jyvec.normalized();
				double jac = Jyvec.norm();

				// evaluate Green's function
				double G = laplace_2d_HSP_kernel()(x, y, nx, ny);
				// multiply by Shape function and Jacobian
				auto N = trial_shape_t::template eval_shape<0>(xi);
				auto F = G * N * jac;
				// evaluate singular part
				auto F0 = (N0 / rho + N1) / rho / twopiJ0;

				// integrate difference numerically
				result.row(i) += (F - F0) * w;
			}

			// add analytic integral of singular part
			result.row(i) += ((1. / (a(0) - xi0(0)) - 1. / (b(0) - xi0(0))) * N0 +
				std::log(std::abs((b(0) - xi0(0)) / (a(0) - xi0(0)))) * N1) / twopiJ0;
		}
		return result;
	}
};

/** \brief Collocational integral of the 2D HSP kernel over a straight line with max. linear Nset */
template <class TestField, class TrialField>
class laplace_2d_HSP_collocation_line
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	typedef typename test_field_t::nset_t test_shape_set_t;
	typedef typename trial_field_t::nset_t trial_shape_set_t;

	typedef typename test_field_t::elem_t elem_t;

	typedef typename test_shape_set_t::xi_t xi_t;

	static size_t const rows = test_shape_set_t::num_nodes;
	static size_t const cols = trial_shape_set_t::num_nodes;

	typedef Eigen::Matrix<double, rows, cols> result_t;

public:
	static result_t eval(elem_t const &elem)
	{
		using namespace boost::math::double_constants;

		double jac = elem.get_normal().norm();

		result_t result;
		result.setZero();

		double a = -1., b = +1.;

		for (size_t i = 0; i < rows; ++i)
		{
			xi_t xi0 = test_shape_set_t::corner_at(i);
			auto N = trial_shape_set_t::template eval_shape<0>(xi0);
			auto dN = trial_shape_set_t::template eval_shape<1>(xi0);
			result.row(i) = N * (-1. / (b - xi0(0)) - -1. / (a - xi0(0))) +
				dN * std::log(std::abs((b - xi0(0)) / (a - xi0(0))));
		}

		return result / (two_pi * jac);
	}
};


/** \brief Collocational singular integral of the 2D Laplace HSP kernel over a straight line with constant shape set */
template <class TestField>
class laplace_2d_HSP_collocation_constant_line
{
	typedef TestField test_field_t;
	typedef typename test_field_t::nset_t test_shape_set_t;
	static size_t const rows = test_shape_set_t::num_nodes;
	typedef Eigen::Matrix<double, rows, 1> result_t;

	typedef line_1_elem elem_t;
	typedef typename elem_t::x_t x_t;

	typedef typename test_shape_set_t::xi_t xi_t;

public:
	/**
	 * \brief Evaluate the integral
	 * \param [in] elem the line element
	 * \param [in] x0 the singular point
	 * \return the integral value
	 */
	static result_t eval(elem_t const &elem)
	{
		using namespace boost::math::double_constants;

		result_t result;
		auto const &C = elem.get_coords();
		for (size_t i = 0; i < rows; ++i)
		{
			x_t x0 = elem.get_x(test_shape_set_t::corner_at(i));
			double d1 = (x0 - C.col(0)).norm(), d2 = (x0 - C.col(1)).norm();
			result(i, 0) = -(1. / d1 + 1. / d2) / two_pi;
		}
		return result;
	}
};



template <class TestField, class TrialField, size_t order>
class laplace_2d_SLP_galerkin_face_general
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	typedef typename test_field_t::nset_t test_shape_t;
	typedef typename trial_field_t::nset_t trial_shape_t;

	static size_t const nTest = test_shape_t::num_nodes;
	static size_t const nTrial = trial_shape_t::num_nodes;

	typedef Eigen::Matrix<double, nTest, nTrial> result_t;

	typedef typename trial_field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;

	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;

	typedef regular_quad_store<domain_t, order> quadrature_t;

public:
	static result_t eval(elem_t const &elem)
	{
		using namespace boost::math::double_constants;

		result_t I = result_t::Zero();

		// loop over u
		for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
		{
			xi_t u = it->get_xi();
			double wu = it->get_w();

			auto Nxi0 = test_shape_t::template eval_shape<0>(u);
			auto Neta0 = trial_shape_t::template eval_shape<0>(u);
			double J0 = elem.get_dx(u).norm();

			// loop over v
			for (auto jt = quadrature_t::quadrature.begin(); jt != quadrature_t::quadrature.end(); ++jt)
			{
				double v = (jt->get_xi()(0) + 1.) / 2.;
				double wv = jt->get_w() / 2.;

				xi_t xi, eta;

				for (int d = 0; d < 2; ++d)
				{
					if (d == 0)
					{
						xi(0) = u(0) + (+1. - u(0)) * v;
						eta(0) = u(0) + (-1. - u(0)) * v;
					}
					else
					{
						xi(0) = u(0) + (-1. - u(0)) * v;
						eta(0) = u(0) + (+1. - u(0)) * v;
					}

					x_t x = elem.get_x(xi);
					x_t y = elem.get_x(eta);
					double G = laplace_2d_SLP_kernel()(x, y);
					double Jxi = elem.get_dx(xi).norm();
					double Jeta = elem.get_dx(xi).norm();
					auto Nxi = test_shape_t::template eval_shape<0>(xi);
					auto Neta = trial_shape_t::template eval_shape<0>(eta);

					// evaluate integrand
					result_t F = G * Jxi * Jeta * (Nxi * Neta.transpose());

					// evaluate integrand's singular part
					double G0 = -std::log(2. * v * J0);
					result_t F0 = G0 * J0 * J0 * (Nxi0 * Neta0.transpose());

					// integrate difference numerically
					I += (F - F0) * (2. * (1 - v)) * wu * wv;
				}
			} // loop over v

			// compensate with analytical inner integral of the regular part
			I += 2. * (1.5 - std::log(2. * J0)) / two_pi * J0 * J0 * (Nxi0 * Neta0.transpose()) * wu;
		} // loop over i

		return I;
	}
};


/** \brief Galerkin face match singular integral of the 2D Laplace SLP kernel over a constant line element */
class laplace_2d_SLP_galerkin_face_constant_line
{
public:
	/**
	* \brief Evaluate the integral
	* \param [in] elem the line element
	* \return the integral value
	*/
	static double eval(line_1_elem const &elem)
	{
		using namespace boost::math::double_constants;

		auto const &C = elem.get_coords();
		double d = (C.col(1) - C.col(0)).norm();	// element length
		return d * d * (1.5 - std::log(d)) / two_pi;
	}
};

/** \brief Galerkin face match singular integral of the 2D Laplace SLP kernel over a linear line element */
class laplace_2d_SLP_galerkin_face_linear_line
{
public:
	/**
	* \brief Evaluate the integral
	* \param [in] elem the line element
	* \param [out] i1 the diagonal elem of the result
	* \param [out] i2 the off-diagonal elem of the result
	*/
	static void eval(line_1_elem const &elem, double &i1, double &i2)
	{
		using namespace boost::math::double_constants;

		auto const &C = elem.get_coords();
		double d = (C.col(1) - C.col(0)).norm();	// element length
		double c = d * d / (8. * pi);
		double lnd = std::log(d);
		i1 = c * (1.75 - lnd);
		i2 = c * (1.25 - lnd);
	}
};

/** \brief Galerkin integral of the 2D SLP kernel over a constant line with edge match */
class laplace_2d_SLP_galerkin_edge_constant_line
{
	static double qfunc(double a, double phi)
	{
		using namespace boost::math::double_constants;

		if (std::abs(phi) < 1e-3)
			return a / (a + 1.);
		double cotphi = std::tan(half_pi - phi);
		return std::atan(a / std::sin(phi) + cotphi) - std::atan(cotphi);
	}

public:
	/** \brief evaluate the integral on two elements */
	static double eval(line_1_elem const &elem1, line_1_elem const &elem2)
	{
		using namespace boost::math::double_constants;

		// get the element corner coordinates
		auto const &C1 = elem1.get_coords();
		auto const &C2 = elem2.get_coords();
		// and side vectors
		auto r1vec = (C1.col(1) - C1.col(0)), r2vec = (C2.col(1) - C2.col(0));
		// and side lengths
		double r1 = r1vec.norm(), r2 = r2vec.norm();
		// and signed angle between them
		double phi = std::asin(r1vec(0) * r2vec(1) - r2vec(0) * r1vec(1)) / (r1 * r2);
		// third side length
		double r3 = std::sqrt(r1 * r1 + 2 * r1 * r2 * std::cos(phi) + r2 * r2);
		return (
			r1 * r2 * (3. - 2. * std::log(r3))
			+ std::cos(phi) * (r1 * r1 * std::log(r1 / r3) + r2 * r2 * std::log(r2 / r3))
			- std::sin(phi) * (r1 * r1 * qfunc(r2 / r1, phi) + r2 * r2 * qfunc(r1 / r2, phi))
			) / (4. * pi);
	}
};

/** \brief Galerkin integral of the 2D DLP kernel over a constant line with edge match */
class laplace_2d_DLP_galerkin_edge_constant_line
{
	static double qfunc(double a, double phi)
	{
		using namespace boost::math::double_constants;

		if (std::abs(phi) < 1e-3)
			return a / (a + 1.);
		double cotphi = std::tan(half_pi - phi);
		return std::atan(a / std::sin(phi) + cotphi) - std::atan(cotphi);
	}

public:
	/** \brief evaluate the integral on two elements */
	static double eval(line_1_elem const &elem1, line_1_elem const &elem2)
	{
		using namespace boost::math::double_constants;

		// get element corners
		auto const &C1 = elem1.get_coords();
		auto const &C2 = elem2.get_coords();
		// get side vectors
		auto r1vec = (C1.col(1) - C1.col(0)), r2vec = (C2.col(1) - C2.col(0));
		// get side lengths
		double r1 = r1vec.norm(), r2 = r2vec.norm();
		// get angle between elements
		double phi = std::asin(r1vec(0) * r2vec(1) - r2vec(0) * r1vec(1)) / (r1 * r2);
		// general expression
		double r3 = std::sqrt(r1 * r1 + 2 * r1 * r2 * std::cos(phi) + r2 * r2);
		double res = (
			r2 * std::cos(phi) * qfunc(r1 / r2, phi)
			- r1 * qfunc(r2 / r1, phi)
			+ r2 * std::sin(phi) * std::log(r2 / r3)
			) / two_pi;

		if (elem1.get_nodes()(0) == elem2.get_nodes()(1))
			res *= -1;

		return res;
	}
};



/** \brief Collocational singular integral of the 3D Laplace SLP kernel over a constant plane element */
class laplace_3d_SLP_collocation_constant_plane
{
public:
	/**
	 * \brief Evaluate the integral
	 * \param [in] elem the line element
	 * \param [in] x0 the singular point
	 * \return the integral value
	 */
	template <class elem_t>
	static double eval(elem_t const &elem, typename elem_t::x_t const &x0)
	{
		using namespace boost::math::double_constants;

		enum { N = elem_t::domain_t::num_corners };
		double r[N], theta[N], alpha[N], result = 0.;
		plane_element_helper(elem, x0, r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result += r[i] * std::sin(alpha[i]) *
			std::log(std::tan((alpha[i] + theta[i]) / 2.) / std::tan(alpha[i] / 2.));

		return result / (4. * pi);
	}
};

/** \brief Collocational singular integral of the 3D Laplace HSP kernel over a constant planar element */
class laplace_3d_HSP_collocation_constant_plane
{
public:
	/**
	 * \brief Evaluate the integral
	 * \param [in] elem the line element
	 * \param [in] x0 the singular point
	 * \return the integral value
	 */
	template <class elem_t>
	static double eval(elem_t const &elem, typename elem_t::x_t const &x0)
	{
		using namespace boost::math::double_constants;

		enum { N = elem_t::domain_t::num_corners };
		double r[N], theta[N], alpha[N], result = 0.;
		plane_element_helper(elem, x0, r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result += (std::cos(alpha[i] + theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));

		return result / (4. * pi);
	}
};


/** \brief collocational singular integral of the 2D SLP kernel over a straight line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value
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
		result = laplace_2d_SLP_collocation_straight_line_second_order<
			TestField, TrialField
		>::eval(trial_field.get_elem());
		return result;
	}
};


/** \brief collocational singular integral of the 2D SLP kernel over a curved line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	!std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] test_field the test field
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_SLP_kernel> const &,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &,
		element_match const &)
	{
		result = laplace_2d_SLP_collocation_general<TestField, TrialField, 10>::eval(
			test_field.get_elem());
		return result;
	}
};


/** \brief Galerkin face-match singular integral of the 2D SLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::match_1d_type,
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
		result(0, 0) = laplace_2d_SLP_galerkin_face_constant_line::eval(trial_field.get_elem());
		return result;
	}
};


/** \brief Galerkin face-match singular integral of the 2D SLP kernel over a linear line
* \tparam TestField the test field type
* \tparam TrialField the trial field type
*/
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
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
		laplace_2d_SLP_galerkin_face_linear_line::eval(result(0, 0), result(0, 1), result(1, 1));
		result(1, 0) = result(0, 1);
		return result;
	}
};

/** \brief Galerkin edge-match singular integral of the 2D SLP kernel over two constant lines
* \tparam TestField the test field type
* \tparam TrialField the trial field type
*/
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::match_0d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value &&
	std::is_same<typename TestField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] test_field the test field
	 * \param [in] trial_field the trial field
	 * \param [in] match the match data
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_SLP_kernel> const &,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &match)
	{
		result(0, 0) = laplace_2d_SLP_galerkin_edge_constant_line::eval(
			test_field.get_elem(), trial_field.get_elem());
		return result;
	}
};


/** \brief Galerkin edge-match singular integral of the 2D DLP kernel over two constant lines
* \tparam TestField the test field type
* \tparam TrialField the trial field type
*/
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_DLP_kernel, TestField, TrialField, match::match_0d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value &&
	std::is_same<typename TestField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] test_field the test field
	 * \param [in] trial_field the trial field
	 * \param [in] match the match data
	 * \return reference to the result matrix
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_DLP_kernel> const &,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &match)
	{
		result(0, 0) = laplace_2d_DLP_galerkin_edge_constant_line::eval(
			test_field.get_elem(), trial_field.get_elem());
		return result;
	}
};


/** \brief collocational singular integral of the 2D HSP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_HSP_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
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
		result = laplace_2d_HSP_collocation_constant_line<TestField>::eval(
			trial_field.get_elem());
		return result;
	}
};


/** \brief collocational singular integral of the 2D HSP kernel over a linear line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_HSP_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	TrialField::nset_t::polynomial_degree == 1
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
		result = laplace_2d_HSP_collocation_line<TestField, TrialField>::eval(
			trial_field.get_elem());
		return result;
	}
};


/** \brief Galerkin face-match singular integral of the 2D HSP kernel over a constant line
* \tparam TestField the test field type
* \tparam TrialField the trial field type
*/
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_HSP_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::general>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
	std::is_same<typename TestField::nset_t, line_0_shape_set>::value &&
	std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_2d_HSP_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &,
		element_match const &)
	{
		return result;
	}
};

/** \brief collocational singular integral of the 3D SLP kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_3d_SLP_kernel, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
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
		result(0, 0) = laplace_3d_SLP_collocation_constant_plane::eval(
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
	laplace_3d_HSP_kernel, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
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
		result(0, 0) = laplace_3d_HSP_collocation_constant_plane::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center());
		return result;
	}
};


/** \brief collocational singular integral of the 3d Gxx kernel over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_3d_Gxx_kernel, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
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
		kernel_base<laplace_3d_Gxx_kernel> const &,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result(0, 0) = -1. * laplace_3d_HSP_collocation_constant_plane::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center());
		return result;
	}
};


/** \brief collocational singular integral of the HSP kernel not over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_3d_HSP_kernel, TestField, TrialField, match::match_2d_type,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	!(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
	>::type
>
{
public:
	/** \brief evaluate singular integral
	 * \tparam result_t the result matrix type
	 * \param [in, out] result reference to the result
	 * \param [in] kernel the kernel instance
	 * \param [in] trial_field the trial and test field
	 */
	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_HSP_kernel> const &kernel,
		field_base<TestField> const &,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		typedef guiggiani<TrialField, laplace_3d_HSP_kernel, 5, 9> guiggiani_t;
		auto const &elem = trial_field.get_elem();
		guiggiani_t gui(elem, kernel.derived());
		for (unsigned r = 0; r < TestField::num_dofs; ++r)
		{
			auto const &xi0 = TestField::nset_t::corner_at(r);
			gui.integrate(result.row(r), xi0, elem.get_normal(xi0));
		}
		return result;
	}
};

} // namespace NiHu

#endif /* LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED */

