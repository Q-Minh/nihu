#ifndef HELMHOLTZ_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED
#define HELMHOLTZ_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "helmholtz_kernel.hpp"
#include "laplace_kernel.hpp"
#include "laplace_nearly_singular_integrals.hpp"
#include "quadrature_store_helper.hpp"

namespace NiHu
{

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
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
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
	helmholtz_3d_DLP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
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
	helmholtz_3d_HSP_kernel<WaveNumber>, TestField, TrialField,
	typename std::enable_if<
		std::is_same<typename get_formalism<TestField, TrialField>::type, formalism::collocational>::value &&
		(std::is_same<typename TrialField::elem_t::lset_t, tria_1_shape_set>::value &&
		std::is_same<typename TrialField::nset_t, tria_0_shape_set>::value)
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

}

#endif
