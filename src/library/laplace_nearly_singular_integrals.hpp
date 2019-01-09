// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2019  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2019  Peter Rucz <rucz@hit.bme.hu>
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

/** \file laplace_nearly_singular_integrals.hpp
 *
 */

#ifndef LAPLACE_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED
#define LAPLACE_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../core/formalism.hpp"
#include "../core/nearly_singular_integral.hpp"
#include "../core/nearly_singular_planar_constant_collocation_shortcut.hpp"
#include "../util/math_functions.hpp"
#include "laplace_kernel.hpp"
#include "plane_element_helper.hpp"
#include "nearly_singular_collocational.hpp"
#include "quadrature_store_helper.hpp"

namespace NiHu
{

/** \brief nearly singular collocational integral of the 3D Laplace SLP kernel over planes */
class laplace_3d_SLP_collocation_constant_plane_nearly_singular
{
	typedef line_domain quadrature_domain_t;
	enum { quadrature_order = 10 };
	typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;

public:
	template <class elem_t>
	static double eval(elem_t const &elem, typename elem_t::x_t const &x0_in)
	{
		typedef typename elem_t::x_t x_t;
		
		enum { N = elem_t::domain_t::num_corners  };

		Eigen::Matrix<double, 3, N> corners;
		for (unsigned i = 0; i < N; ++i)
			corners.col(i) = elem.get_x(elem_t::domain_t::get_corner(i));

		auto T = plane_elem_transform(corners.col(1)-corners.col(0), corners.col(2)-corners.col(0));
		auto Tdec = T.partialPivLu();
		corners = Tdec.solve(corners);
		x_t x0 = Tdec.solve(x0_in);
		
		double theta_lim[N];
		double ref_distance[N];
		double theta0[N];
		
		plane_elem_helper_mid(corners, x0, ref_distance, theta_lim, theta0);
		
		
		double z = x0(2) - corners(2, 0);
		
		double result = 0.0;

		// loop over triange sides
		for (unsigned i = 0; i < N; ++i)
		{
			double th1 = theta_lim[i];
			double th2 = theta_lim[(i+1)%N];

			// angle check
			if (std::abs(th2 - th1) > M_PI)
			{
				if (th2 > th1) th1 += 2 * M_PI;
				else th2 += 2*M_PI;
			}
				
			double jac = (th2 - th1)/2.;
			
			// loop over quadrature points
			for (auto it = quadrature_t::quadrature.begin();
				it != quadrature_t::quadrature.end(); ++it)
			{
				double xi = it->get_xi()(0);
				double w = it->get_w();
				
				// apply shape function to get angle
				double theta = th1 * (1.-xi)/2. + th2 * (1.+xi)/2.;
				
				double R = ref_distance[i] / std::cos(theta - theta0[i]);
				double r = std::sqrt(R*R + z*z);
				
				// although z is constant, it can not be extracted from the 
				// numerical integration, as int dtheta is 2pi or 0
				result += (r - std::abs(z)) * jac * w;
			}
		}
		
		return result / (4. * M_PI);
	}
};


class laplace_3d_DLP_collocation_constant_plane_nearly_singular
{
	typedef line_domain quadrature_domain_t;
	enum { quadrature_order = 10 };
	typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;

public:
	template <class elem_t>
	static double eval(elem_t const &elem, typename elem_t::x_t const &x0_in)
	{
		typedef typename elem_t::x_t x_t;
		
		enum { N = elem_t::domain_t::num_corners  };
		auto corners = elem.get_coords();
		
		auto T = plane_elem_transform(corners.col(1)-corners.col(0), corners.col(2)-corners.col(0));
		auto Tdec = T.partialPivLu();
		corners = Tdec.solve(corners);
		x_t x0 = Tdec.solve(x0_in);
		double z = x0(2) - corners(2,0);
		
		double theta_lim[N];
		double ref_distance[N];
		double theta0[N];
		
		plane_elem_helper_mid(corners, x0, ref_distance, theta_lim, theta0);
		
		double result = 0.0;

		// loop over triange sides
		for (unsigned i = 0; i < N; ++i)
		{
			double th1 = theta_lim[i];
			double th2 = theta_lim[(i+1)%N];

			// angle check
			if (std::abs(th2 - th1) > M_PI)
			{
				if (th2 > th1) th1 += 2 * M_PI;
				else th2 += 2*M_PI;
			}
				
			double jac = (th2 - th1)/2.;
			
			// loop over quadrature points
			for (auto it = quadrature_t::quadrature.begin();
				it != quadrature_t::quadrature.end(); ++it)
			{
				double xi = it->get_xi()(0);
				double w = it->get_w();
				
				// apply shape function to get angle
				double theta = th1 * (1.-xi)/2. + th2 * (1.+xi)/2.;
				
				double R = ref_distance[i] / std::cos(theta - theta0[i]);
				double r = std::sqrt(R*R + z*z);
				
				// although z is constant, it can not be extracted from the 
				// numerical integration, as int dtheta is 2pi or 0
				result += (1 - std::abs(z) / r) * jac * w;
			}
		}
		
		return sgn(z) * result / (4. * M_PI);
	}
};


class laplace_3d_DLPt_collocation_constant_plane_nearly_singular
{
	typedef line_domain quadrature_domain_t;
	enum { quadrature_order = 10 };
	typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;

public:
	template <class elem_t>
	static double eval(
		elem_t const &elem,
		typename elem_t::x_t const &x0_in,
		typename elem_t::x_t const &nx_in)
	{
		typedef typename elem_t::x_t x_t;

		enum { N = elem_t::domain_t::num_corners };
		auto corners = elem.get_coords();

		auto T = plane_elem_transform(corners.col(1) - corners.col(0),
			corners.col(2) - corners.col(0));
		auto Tdec = T.partialPivLu();
		corners = Tdec.solve(corners);
		x_t x0 = Tdec.solve(x0_in);
		x_t nx = Tdec.solve(nx_in);

		double z = x0(2) - corners(2,0);
		
		double theta_lim[N];
		double ref_distance[N];
		double theta0[N];
		
		plane_elem_helper_mid(corners, x0, ref_distance, theta_lim, theta0);
		
		double result = 0.0;

		// loop over triange sides
		for (unsigned i = 0; i < N; ++i)
		{
			double th1 = theta_lim[i];
			double th2 = theta_lim[(i+1)%N];

			// angle check
			if (std::abs(th2 - th1) > M_PI)
			{
				if (th2 > th1) th1 += 2 * M_PI;
				else th2 += 2*M_PI;
			}
				
			double jac = (th2 - th1)/2.;
			
			// loop over quadrature points
			for (auto it = quadrature_t::quadrature.begin();
				it != quadrature_t::quadrature.end(); ++it)
			{
				double xi = it->get_xi()(0);
				double w = it->get_w();
				
				// apply shape function to get angle
				double theta = th1 * (1.-xi)/2. + th2 * (1.+xi)/2.;
				
				double R = ref_distance[i] / std::cos(theta - theta0[i]);
				double r = std::sqrt(R*R + z*z);
				
				x_t rvec;
				rvec(0) = R * std::cos(theta);
				rvec(1) = R * std::sin(theta);
				rvec(2) = -z;
				
				// although z is constant, it can not be extracted from the 
				// numerical integration, as int dtheta is 2pi or 0
				result += (
					-nx.dot(rvec) / r +
					(nx(0) * std::cos(theta) + nx(1) * std::sin(theta)) * std::log((R + r))
					- nx(2) * sgn(z)
				) * jac * w;
			}
		}

		return result / (4.*M_PI);
	}
};


class laplace_3d_HSP_collocation_constant_plane_nearly_singular
{
	typedef line_domain quadrature_domain_t;
	enum { quadrature_order = 40 };
	typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;

public:
	template <class elem_t>
	static double eval(
		elem_t const &elem, 
		typename elem_t::x_t const &x0_in, 
		typename elem_t::x_t const &nx_in)
	{
		typedef typename elem_t::x_t x_t;
		
		enum { N = elem_t::domain_t::num_corners  };
		auto corners = elem.get_coords();
		
		auto T = plane_elem_transform(corners.col(1)-corners.col(0), corners.col(2)-corners.col(0));
		auto Tdec = T.partialPivLu();
		corners = Tdec.solve(corners);
		x_t x0 = Tdec.solve(x0_in);
		x_t nx = Tdec.solve(nx_in);
		
		double theta_lim[N];
		double ref_distance[N];
		double theta0[N];
	
		
		plane_elem_helper_mid(corners, x0, ref_distance, theta_lim, theta0);
		
		double result = 0.0;
		double z = x0(2) - corners(2,0);

		// loop over triange sides
		for (unsigned i = 0; i < N; ++i)
		{
			double th1 = theta_lim[i];
			double th2 = theta_lim[(i+1)%N];

			// angle check
			if (std::abs(th2 - th1) > M_PI)
			{
				if (th2 > th1) th1 += 2 * M_PI;
				else th2 += 2*M_PI;
			}
				
			double jac = (th2 - th1)/2.;
			
			// loop over quadrature points
			for (auto it = quadrature_t::quadrature.begin();
				it != quadrature_t::quadrature.end(); ++it)
			{
				double xi = it->get_xi()(0);
				double w = it->get_w();
				
				// apply shape function to get angle
				double theta = th1 * (1.-xi)/2. + th2 * (1.+xi)/2.;
				
				double R = ref_distance[i] / std::cos(theta - theta0[i]);
				double r = std::sqrt(R*R + z*z);
				
				x_t rvec;
				rvec(0) = R * std::cos(theta);
				rvec(1) = R * std::sin(theta);
				rvec(2) = -z;
				
				double rdnx = -rvec.dot(nx) / r;
				double integrand;
				
				/** \todo Check z->0 limiting case */
				if (std::abs(z) > 1e-12)
					integrand = -R*R / (r*r) / z * rdnx;
				else
					integrand = -nx(2) / R;
				
				result += integrand * jac * w;
			}
		}

		return result / (4.*M_PI);
	}
};

/**
 * \brief Class enabling the specialisation for 3D SLP Laplace kernel
 * \tparam TestField Test field type
 * \tparam TrialField trial field type
 * \details
 * The specialisation is enabled for the collocational formalism if the element 
 * is a linear triangle and the field is constant.
 */
template <class TestField, class TrialField>
class nearly_singular_integral<
	laplace_3d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
	>::type
>
{
public:
	/**
	 * \brief Check if singular evaluation is needed 
	 * \todo needed function should be implemented in a base class
	 */
	static bool needed(
		kernel_base<laplace_3d_SLP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		double const limit = 1.5;
		
		// distance between element centres
		double d = (test_field.get_elem().get_center() - trial_field.get_elem().get_center()).norm();
		double R = trial_field.get_elem().get_linear_size_estimate();
		return d/R < limit;
	}

	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_SLP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		result(0) = laplace_3d_SLP_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center());
		return result;
	}
};

template <class TestField, class TrialField>
class nearly_singular_integral<
	laplace_3d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	! (
	(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value) )
	>::type
>
{
public:
	static bool needed(
		kernel_base<laplace_3d_SLP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
	)
	{
		double const limit = 1.5;

		// distance between element centres
		double d = (test_field.get_elem().get_center() - trial_field.get_elem().get_center()).norm();
		double R = trial_field.get_elem().get_linear_size_estimate();
		return d / R < limit;
	}

	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_SLP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
	)
	{
		nearly_singular_collocational<TrialField, laplace_3d_SLP_kernel, 5, 5> nsg(trial_field, kernel);
		laplace_3d_SLP_kernel::test_input_t tsi(test_field.get_elem(), TestField::elem_t::domain_t::get_center());
		nsg.integrate(result, tsi);
		return result;
	}
};




template <class TestField, class TrialField>
class nearly_singular_integral<
	laplace_3d_DLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
	>::type
>
{
public:
	static bool needed(
		kernel_base<laplace_3d_DLP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		double const limit = 1.5;
		
		// distance between element centres
		double d = (test_field.get_elem().get_center() - trial_field.get_elem().get_center()).norm();
		double R = trial_field.get_elem().get_linear_size_estimate();
		return d/R < limit;
	}

	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_DLP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		result(0) = laplace_3d_DLP_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center());
		return result;
	}
};


template <class TestField, class TrialField>
class nearly_singular_integral<
	laplace_3d_DLPt_kernel, TestField, TrialField,
	typename std::enable_if<
	std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
	(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
	>::type
>
{
public:
	static bool needed(
		kernel_base<laplace_3d_DLPt_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
	)
	{
		double const limit = 1.5;

		// distance between element centres
		double d = (test_field.get_elem().get_center() - trial_field.get_elem().get_center()).norm();
		double R = trial_field.get_elem().get_linear_size_estimate();
		return d / R < limit;
	}

	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_DLPt_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
	)
	{
		result(0) = laplace_3d_DLPt_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center(),
			test_field.get_elem().get_normal().normalized());
		return result;
	}
};



template <class TestField, class TrialField>
class nearly_singular_integral<
	laplace_3d_HSP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
	>::type
>
{
public:
	static bool needed(
		kernel_base<laplace_3d_HSP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		double const limit = 1.5;
		
		// distance between element centres
		double d = (test_field.get_elem().get_center() - trial_field.get_elem().get_center()).norm();
		double R = trial_field.get_elem().get_linear_size_estimate();
		return d/R < limit;
	}

	template <class result_t>
	static result_t &eval(
		result_t &result,
		kernel_base<laplace_3d_HSP_kernel> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		result(0) = laplace_3d_HSP_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center(),
			test_field.get_elem().get_normal().normalized());
		return result;
	}
};

template <class Elem>
class nearly_singular_planar_constant_collocation_shortcut<laplace_3d_SLP_kernel, Elem>
{
public:
	typedef Elem elem_t;
	typedef laplace_3d_SLP_kernel::result_t res_t;

	static res_t eval(
		laplace_3d_SLP_kernel::test_input_t const &test_input,
		elem_t const &elem)
	{
		return laplace_3d_SLP_collocation_constant_plane_nearly_singular::eval(
			elem, test_input.get_x());
	}
};


} // end of namespace NiHu

#endif
