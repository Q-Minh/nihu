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

/** \file laplace_singular_integrals.hpp
 * \brief Analytical expressions for the singular integrals of Laplace kernels
 * \details analytical expression for the Laplace kernels over plane elements
 */
#ifndef LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED
#define LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../core/match_types.hpp"
#include "../core/singular_integral_shortcut.hpp"
#include "../core/quadrature.hpp"
#include "lib_element.hpp"
#include "laplace_kernel.hpp"
#include "normal_derivative_singular_integrals.hpp"
#include "plane_element_helper.hpp"
#include "guiggiani_1992.hpp"
#include "lib_element.hpp"

#include <iostream>

namespace NiHu
{
	
/** \brief store-wrapper of a statically stored quadrature */
template <class domain_t, size_t order>
struct regular_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<domain_t> const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <class domain_t, size_t order>
gaussian_quadrature<domain_t> const regular_quad_store<domain_t, order>::quadrature(order);

class log_gauss_quadrature
	: public std::vector<quadrature_elem<Eigen::Matrix<double, 1, 1>, double> >
{
	typedef Eigen::Matrix<double, 1, 1> xi_t;
	typedef double w_t;
	typedef quadrature_elem<xi_t, w_t> quadrature_elem_t;
	typedef std::vector<quadrature_elem_t> base_t;
	
public:
	log_gauss_quadrature(size_t points)
	{
		base_t::reserve(points);
		
		switch (points)
		{
		case 1:
			push_back(quadrature_elem_t(xi_t::Constant(0.25), 1.));
			break;
		case 2:
			push_back(quadrature_elem_t(xi_t::Constant(0.602276908118738), 0.281460680969616));
			push_back(quadrature_elem_t(xi_t::Constant(0.112008806166976), 0.718539319030385));
			break;
		case 3:
			push_back(quadrature_elem_t(xi_t::Constant(0.766880303938942), 0.094615406566148));
			push_back(quadrature_elem_t(xi_t::Constant(0.368997063715618), 0.391980041201488));
			push_back(quadrature_elem_t(xi_t::Constant(0.063890793087325), 0.513404552232364));
			break;
		default:
			throw std::out_of_range("unsupported log gauss degree");
			break;
		}
	}
};


template <unsigned points>
struct singular_quad_store
{
	static log_gauss_quadrature const quadrature;
};

template <unsigned points>
log_gauss_quadrature const singular_quad_store<points>::quadrature(points);
	

template <class TestField, class TrialField>
class laplace_2d_SLP_collocation_line
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	
	typedef typename test_field_t::nset_t test_shape_set_t;
	typedef typename trial_field_t::nset_t trial_shape_set_t;
	
	typedef typename test_field_t::elem_t elem_t;
	
	typedef typename test_shape_set_t::xi_t xi_t;
	
	static size_t const rows = test_shape_set_t::num_nodes;
	static size_t const cols = trial_shape_set_t::num_nodes;
	static size_t const order = trial_shape_set_t::polynomial_order;
	
	typedef Eigen::Matrix<double, rows, cols> result_t;
	
	typedef regular_quad_store<typename elem_t::domain_t, order> reg_t;
	typedef singular_quad_store<(order + 1 + 1)/2> sing_t;
		
public:
	static result_t eval(elem_t const &elem)
	{
		// get Jacobian
		double jac = elem.get_normal().norm();
		
		result_t result;
		result.setZero();
		
		
		// traverse collocational points
		for (size_t i = 0; i < rows; ++i)
		{
			xi_t xi0 = test_shape_set_t::corner_at(i);
			double rho1 = (1. + xi0(0));
			double rho2 = (1. - xi0(0));
			
			for (auto it = reg_t::quadrature.begin(); it != reg_t::quadrature.end(); ++it)
			{
				// get quadrature location and weight
				xi_t const &zeta = it->get_xi();
				double w = it->get_w();
				
				xi_t eta = (xi_t::Constant(1.0) + zeta)/2.;
				
				result.row(i) +=
				(trial_shape_set_t::template eval_shape<0>(xi0 - rho1*eta) * rho1/2. * std::log(1./rho1) +
				trial_shape_set_t::template eval_shape<0>(xi0 + rho2*eta) * rho2/2. * std::log(1./rho2) +
				trial_shape_set_t::template eval_shape<0>(zeta) * std::log(1./jac)).transpose() * w;
			}
			
			for (auto it = sing_t::quadrature.begin(); it != sing_t::quadrature.end(); ++it)
			{
				// get quadrature location and weight
				xi_t const &eta = it->get_xi();
				double w = it->get_w();
				
				result.row(i) +=
				(trial_shape_set_t::template eval_shape<0>(xi0 - rho1*eta) * rho1 +
				trial_shape_set_t::template eval_shape<0>(xi0 + rho2*eta) * rho2).transpose() * w;
			}
		}
		
		return result * jac / (2. *M_PI);
	}
};


/** \brief Collocational singular integral of the 2D Laplace SLP kernel over a constant line element */
class laplace_2d_SLP_collocation_constant_line
{
public:
    /**
     * \brief Evaluate the integral
     * \param [in] elem the line element
     * \param [in] x0 the singular point
     * \return the integral value
     */
	static double eval(line_1_elem const &elem, line_1_elem::x_t const &x0)
	{
		auto const &C = elem.get_coords();
		double d1 = (x0 - C.col(0)).norm(), d2 = (x0 - C.col(1)).norm();
		return (d1 * (1. - std::log(d1)) + d2 * (1. - std::log(d2))) / (2.*M_PI);
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
		auto const &C = elem.get_coords();
		double d = (C.col(1) - C.col(0)).norm();	// element length
		return d*d*(1.5-std::log(d)) / (2.*M_PI);
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
		auto const &C = elem.get_coords();
		double d = (C.col(1) - C.col(0)).norm();	// element length
		double c = d*d/(8.*M_PI);
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
    	if (std::abs(phi) < 1e-3)
    		return a / (a + 1.);
    	double cotphi = std::tan(M_PI/2.-phi);
        return std::atan(a/std::sin(phi) + cotphi) - std::atan(cotphi);
    }

public:
    /** \brief evaluate the integral on two elements */
    static double eval(line_1_elem const &elem1, line_1_elem const &elem2)
    {
		// get the element corner coordinates
        auto const &C1 = elem1.get_coords();
        auto const &C2 = elem2.get_coords();
        // and side vectors
        auto r1vec = (C1.col(1) - C1.col(0)), r2vec = (C2.col(1) - C2.col(0));
        // and side lengths
        double r1 = r1vec.norm(), r2 = r2vec.norm();
        // and signed angle between them
        double phi = std::asin(r1vec(0)*r2vec(1)-r2vec(0)*r1vec(1))/(r1*r2);
		// third side length
		double r3 = std::sqrt(r1*r1 + 2*r1*r2*std::cos(phi) + r2*r2);
		return (
			r1*r2*(3.-2.*std::log(r3))
			+ std::cos(phi) * (r1*r1*std::log(r1/r3)+r2*r2*std::log(r2/r3))
			- std::sin(phi) * (r1*r1*qfunc(r2/r1, phi) + r2*r2*qfunc(r1/r2, phi))
		) / (4.*M_PI);
    }
};

/** \brief Galerkin integral of the 2D DLP kernel over a constant line with edge match */
class laplace_2d_DLP_galerkin_edge_constant_line
{
    static double qfunc(double a, double phi)
    {
    	if (std::abs(phi) < 1e-3)
    		return a / (a + 1.);
    	double cotphi = std::tan(M_PI/2.-phi);
        return std::atan(a/std::sin(phi) + cotphi) - std::atan(cotphi);
    }

public:
    /** \brief evaluate the integral on two elements */
    static double eval(line_1_elem const &elem1, line_1_elem const &elem2)
    {
		// get element corners
        auto const &C1 = elem1.get_coords();
        auto const &C2 = elem2.get_coords();
        // get side vectors
        auto r1vec = (C1.col(1) - C1.col(0)), r2vec = (C2.col(1) - C2.col(0));
        // get side lengths
        double r1 = r1vec.norm(), r2 = r2vec.norm();
        // get angle between elements
        double phi = std::asin(r1vec(0)*r2vec(1)-r2vec(0)*r1vec(1))/(r1*r2);
		// general expression
		double r3 = std::sqrt(r1*r1 + 2*r1*r2*std::cos(phi) + r2*r2);
		double res = (
			r2*std::cos(phi) * qfunc(r1/r2, phi)
			- r1 * qfunc(r2/r1, phi)
			+ r2*std::sin(phi) * std::log(r2/r3)
		) / (2.*M_PI);
		
		if (elem1.get_nodes()(0) == elem2.get_nodes()(1))
			res *= -1;
		
		return res;
    }
};


/** \brief Collocational integral of the 2D HSP kernel over a line with max. linear Nset */
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
		double jac = elem.get_normal().norm();
		
		result_t result;
		result.setZero();

		double a = -1., b = +1.;
		
		for (size_t i = 0; i < rows; ++i)
		{
			xi_t xi0 = test_shape_set_t::corner_at(i);
			auto N = trial_shape_set_t::template eval_shape<0>(xi0);
			auto dN = trial_shape_set_t::template eval_shape<1>(xi0);
			result.row(i) = N * (-1./(b-xi0(0)) - -1./(a-xi0(0)) ) +
				dN * std::log( std::abs( (b-xi0(0))/(a-xi0(0)) ) );
		}
		
		return result / (2. *M_PI * jac);
	}
};


/** \brief Collocational singular integral of the 2D Laplace HSP kernel over a constant line element */
class laplace_2d_HSP_collocation_constant_line
{
public:
    /**
     * \brief Evaluate the integral
     * \param [in] elem the line element
     * \param [in] x0 the singular point
     * \return the integral value
     */
	static double eval(line_1_elem const &elem, line_1_elem::x_t const &x0)
	{
		auto const &C = elem.get_coords();
		double d1 = (x0 - C.col(0)).norm(), d2 = (x0 - C.col(1)).norm();
		return -(1. / d1 + 1. / d2) / (2.*M_PI);
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
		enum { N = elem_t::domain_t::num_corners  };
		double r[N], theta[N], alpha[N], result = 0.;
		plane_element_helper(elem, x0, r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result += r[i] * std::sin(alpha[i]) *
			std::log(std::tan((alpha[i] + theta[i]) / 2.) / std::tan(alpha[i] / 2.));

		return result / (4.*M_PI);
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
		enum { N = elem_t::domain_t::num_corners };
		double r[N], theta[N], alpha[N], result = 0.;
		plane_element_helper(elem, x0, r, theta, alpha);

		for (unsigned i = 0; i < N; ++i)
			result += (std::cos(alpha[i] + theta[i]) - std::cos(alpha[i])) / (r[i] * std::sin(alpha[i]));

		return result / (4.*M_PI);
	}
};


/** \brief collocational singular integral of the 2D SLP kernel over a constant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::match_1d_type,
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


/** \brief collocational singular integral of the 2D SLP kernel over a nonconstant line
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class singular_integral_shortcut<
	laplace_2d_SLP_kernel, TestField, TrialField, match::match_1d_type,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		std::is_same<typename TrialField::elem_t::lset_t, line_1_shape_set>::value &&
		!std::is_same<typename TrialField::nset_t, line_0_shape_set>::value
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
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field,
		element_match const &)
	{
		result = laplace_2d_SLP_collocation_line<TestField, TrialField>::eval(
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
		result(0,0) = laplace_2d_SLP_galerkin_edge_constant_line::eval(
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
		result(0,0) = laplace_2d_DLP_galerkin_edge_constant_line::eval(
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
		result(0,0) = laplace_2d_HSP_collocation_constant_line::eval(
			trial_field.get_elem(),
			trial_field.get_elem().get_center());
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

#endif // LAPLACE_SINGULAR_INTEGRALS_HPP_INCLUDED

