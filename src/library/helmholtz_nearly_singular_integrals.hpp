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

/** 
 * \file helmholtz_nearly_singular_integrals.hpp
 * \brief Nearly singular integrals for Helmholtz kernels
 */

#ifndef HELMHOLTZ_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED
#define HELMHOLTZ_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "helmholtz_kernel.hpp"
#include "laplace_kernel.hpp"
#include "laplace_nearly_singular_integrals.hpp"
#include "nearly_singular_collocational.hpp"
#include "quadrature_store_helper.hpp"

#include "../core/nearly_singular_planar_constant_collocation_shortcut.hpp"

namespace NiHu
{

template <class TestField, class TrialField>
struct is_constant_tria : std::integral_constant<
	bool,
	std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value
	&&
	std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value
> {};

/** 
 * \brief Nearly singular collocational integral of the 3D Helmholtz SLP kernel over planes
 * \details
 * The singular static part is redirected to the corresponding Laplace kernel-
 * The regular dynamic part is integrated numerically using a high-order regular
 * quadrature.
 */
class helmholtz_3d_SLP_collocation_constant_plane_nearly_singular
{
public:
	template <class elem_t, class wavenumber_t>
	static std::complex<double> eval(
		elem_t const &elem,
		typename elem_t::x_t const &x0,
		wavenumber_t const &k)
	{
		typedef typename elem_t::domain_t quadrature_domain_t;
		enum { quadrature_order = 14 };
		typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;
		typedef typename elem_t::x_t x_t;
		
		helmholtz_3d_SLP_kernel<wavenumber_t> ktotal(k);
		laplace_3d_SLP_kernel kstatic;
		
		double res_static = laplace_3d_SLP_collocation_constant_plane_nearly_singular::eval(elem, x0);
		
		std::complex<double> res_dynamic = 0.0;
		
		for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
		{
			auto const &xi = it->get_xi();
			double w = it->get_w();
			
			x_t y = elem.get_x(xi);
			double jac = elem.get_normal(xi).norm();
			
			res_dynamic += (
				ktotal(x0, y) - kstatic(x0, y)
			) * w * jac;
		}

		return res_static + res_dynamic;
	}
};


/** \brief nearly singular collocational integral of the 3D Helmholtz DLP kernel over planes
 * The singular static part is redirected to the corresponding Laplace kernel
 * The regular dynamic part is integrated numerically using a regular quadrature
 */
class helmholtz_3d_DLP_collocation_constant_plane_nearly_singular
{
public:
	template <class elem_t, class wavenumber_t>
	static std::complex<double> eval(
		elem_t const &elem,
		typename elem_t::x_t const &x0,
		wavenumber_t const &k)
	{
		typedef typename elem_t::domain_t quadrature_domain_t;
		enum { quadrature_order = 14 };
		typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;
		typedef typename elem_t::x_t x_t;
		
		helmholtz_3d_DLP_kernel<wavenumber_t> ktotal(k);
		laplace_3d_DLP_kernel kstatic;
		
		double res_static = laplace_3d_DLP_collocation_constant_plane_nearly_singular::eval(elem, x0);
		
		std::complex<double> res_dynamic = 0.0;
		
		x_t ny = elem.get_normal().normalized();
		
		for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
		{
			auto const &xi = it->get_xi();
			double w = it->get_w();
			
			x_t y = elem.get_x(xi);
			double jac = elem.get_normal(xi).norm();
			
			res_dynamic += (
				ktotal(x0, y, ny) - kstatic(x0, y, ny)
			) * w * jac;
		}

		return res_static + res_dynamic;
	}
};

/** 
 * \brief nearly singular collocational integral of the 3D Helmholtz DLPt kernel over planes
 * \details 
 * The singular static part is redirected to the corresponding Laplace kernel
 * The regular dynamic part is integrated numerically using a regular quadrature
 */
class helmholtz_3d_DLPt_collocation_constant_plane_nearly_singular
{
public:
	template <class elem_t, class wavenumber_t>
	static std::complex<double> eval(
		elem_t const &elem,
		typename elem_t::x_t const &x0,
		typename elem_t::x_t const &nx,
		wavenumber_t const &k)
	{
		typedef typename elem_t::domain_t quadrature_domain_t;
		enum { quadrature_order = 14 };
		typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;
		typedef typename elem_t::x_t x_t;
		
		helmholtz_3d_DLPt_kernel<wavenumber_t> ktotal(k);
		laplace_3d_DLPt_kernel kstatic;
		
		double res_static = laplace_3d_DLPt_collocation_constant_plane_nearly_singular::eval(elem, x0, nx);
		
		std::complex<double> res_dynamic = 0.0;
		
		for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
		{
			auto const &xi = it->get_xi();
			double w = it->get_w();
			
			x_t y = elem.get_x(xi);
			double jac = elem.get_normal(xi).norm();
			
			res_dynamic += (
				ktotal(x0, y, nx) - kstatic(x0, y, nx)
			) * w * jac;
		}

		return res_static + res_dynamic;
	}
};


/** \brief nearly singular collocational integral of the 3D Helmholtz HSP kernel over planes
 * The singular static part is redirected to the corresponding Laplace kernel
 * The regular dynamic part is integrated numerically using a regular quadrature
 */
class helmholtz_3d_HSP_collocation_constant_plane_nearly_singular
{
public:
	template <class elem_t, class wavenumber_t>
	static std::complex<double> eval(
		elem_t const &elem,
		typename elem_t::x_t const &x0,
		typename elem_t::x_t const &nx,
		wavenumber_t const &k)
	{
		typedef typename elem_t::domain_t quadrature_domain_t;
		enum { quadrature_order = 14 };
		typedef regular_quad_store<quadrature_domain_t, quadrature_order> quadrature_t;
		typedef typename elem_t::x_t x_t;
		
		helmholtz_3d_HSP_kernel<wavenumber_t> ktotal(k);
		laplace_3d_HSP_kernel kstatic;
		
		double res_static = laplace_3d_HSP_collocation_constant_plane_nearly_singular::eval(elem, x0, nx);
		
		std::complex<double> res_dynamic = 0.0;
		x_t ny = elem.get_normal().normalized();
		
		for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
		{
			auto const &xi = it->get_xi();
			double w = it->get_w();
			
			x_t y = elem.get_x(xi);
			double jac = elem.get_normal(xi).norm();
			
			res_dynamic += (
				ktotal(x0, y, nx, ny) - kstatic(x0, y, nx, ny)
			) * w * jac;
		}

		return res_static + res_dynamic;
	}
};



template <class TestField, class TrialField, class WaveNumber>
class nearly_singular_integral<
	helmholtz_3d_SLP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		is_collocational<TestField, TrialField>::value 
		&&
		is_constant_tria<TestField, TrialField>::value
	>::type
>
{
	typedef helmholtz_3d_SLP_kernel<WaveNumber> kernel_t;
public:
	static bool needed(
		kernel_base<kernel_t> const &kernel,
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
		kernel_base<kernel_t> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		result(0) = helmholtz_3d_SLP_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center(),
			kernel.derived().get_wave_number());
		return result;
	}
};


template <class TestField, class TrialField, class WaveNumber>
class nearly_singular_integral<
	helmholtz_3d_SLP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		is_collocational<TestField, TrialField>::value 
		&&
		!(is_constant_tria<TestField, TrialField>::value)
	>::type
>
{
	typedef helmholtz_3d_SLP_kernel<WaveNumber> kernel_t;
	typedef typename kernel_t::test_input_t test_input_t;
	typedef TestField test_field_t;
	typedef typename test_field_t::nset_t test_nset_t;
	static unsigned const num_test_nodes = test_nset_t::num_nodes;
	
public:
	static bool needed(
		kernel_base<kernel_t> const &kernel,
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
		kernel_base<kernel_t> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		typedef nearly_singular_collocational<TrialField, kernel_t, 15, 15> nsc_t;
		nsc_t nsc(trial_field, kernel);
		
		for (unsigned i = 0; i < num_test_nodes; ++i)
		{
			test_input_t tsti(test_field.get_elem(), test_nset_t::corner_at(i));
			nsc.integrate(result.row(i), tsti);
		}
		
		return result;
	}
};


template <class TestField, class TrialField, class WaveNumber>
class nearly_singular_integral<
	helmholtz_3d_DLP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		is_collocational<TestField, TrialField>::value 
		&&
		is_constant_tria<TestField, TrialField>::value
	>::type
>
{
	typedef helmholtz_3d_DLP_kernel<WaveNumber> kernel_t;
public:
	static bool needed(
		kernel_base<kernel_t> const &kernel,
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
		kernel_base<kernel_t> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		result(0) = helmholtz_3d_DLP_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center(),
			kernel.derived().get_wave_number());
		return result;
	}
};

template <class TestField, class TrialField, class WaveNumber>
class nearly_singular_integral<
	helmholtz_3d_DLPt_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		is_collocational<TestField, TrialField>::value 
		&&
		is_constant_tria<TestField, TrialField>::value
	>::type
>
{
	typedef helmholtz_3d_DLPt_kernel<WaveNumber> kernel_t;
public:
	static bool needed(
		kernel_base<kernel_t> const &kernel,
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
		kernel_base<kernel_t> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		result(0) = helmholtz_3d_DLPt_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center(),
			test_field.get_elem().get_normal().normalized(),
			kernel.derived().get_wave_number());
		return result;
	}
};


template <class TestField, class TrialField, class WaveNumber>
class nearly_singular_integral<
	helmholtz_3d_HSP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		is_collocational<TestField, TrialField>::value 
		&&
		is_constant_tria<TestField, TrialField>::value
	>::type
>
{
	typedef helmholtz_3d_HSP_kernel<WaveNumber> kernel_t;
public:
	static bool needed(
		kernel_base<kernel_t> const &kernel,
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
		kernel_base<kernel_t> const &kernel,
		field_base<TestField> const &test_field,
		field_base<TrialField> const &trial_field
		)
	{
		result(0) = helmholtz_3d_HSP_collocation_constant_plane_nearly_singular::eval(
			trial_field.get_elem(),
			test_field.get_elem().get_center(),
			test_field.get_elem().get_normal().normalized(),
			kernel.derived().get_wave_number());
		return result;
	}
};


template <class Elem, class WaveNumber>
class nearly_singular_planar_constant_collocation_shortcut<helmholtz_3d_SLP_kernel<WaveNumber>, Elem>
{
public:
	typedef Elem elem_t;
	typedef helmholtz_3d_SLP_kernel<WaveNumber> kernel_t;
	typedef typename kernel_t::test_input_t test_input_t;
	typedef typename kernel_t::result_t res_t;

	static res_t eval(
		test_input_t const &test_input,
		elem_t const &elem,
		kernel_base<kernel_t> const &kernel)
	{
		return helmholtz_3d_SLP_collocation_constant_plane_nearly_singular::eval(
			elem,
			test_input.get_x(),
			kernel.derived().get_wave_number());
	}
};





} // end of namespace NiHu

#endif /* HELMHOLTZ_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED */
