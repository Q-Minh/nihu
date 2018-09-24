#ifndef LAPLACE_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED
#define LAPLACE_NEARLY_SINGULAR_INTEGRALS_HPP_INCLUDED

#include "../core/formalism.hpp"
#include "../core/nearly_singular_integral.hpp"
#include "plane_element_helper.hpp"
#include "quadrature_store_helper.hpp"

namespace NiHu
{

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
		auto corners = elem.get_coords();
		
		auto T = plane_elem_transform(corners.col(1)-corners.col(0), corners.col(2)-corners.col(0));
		auto Tdec = T.partialPivLu();
		corners = Tdec.solve(corners);
		x_t x0 = Tdec.solve(x0_in);
		
		double result = 0;

		for (unsigned i = 0; i < N; ++i)
		{
			auto const &c1 = corners.col(i);
			auto const &c2 = corners.col((i+1) % N);
			
			x_t dyxi = (c2-c1)/2.;
			
			for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
			{
				double xi = it->get_xi()(0);
				double w = it->get_w();
				
				x_t y = c1 * (1.-xi)/2. + c2 * (1.+xi)/2.;
				
				x_t rvec = y - x0;
				double r = rvec.norm();
				double z = -rvec(2);
				// square of lateral radius
				double R2 = rvec(0)*rvec(0) + rvec(1)*rvec(1);
				
				double integrand = r-std::abs(z);
				
				double dtheta = (rvec(0) * dyxi(1) - rvec(1) * dyxi(0)) / R2;
				
				result += integrand * dtheta * w;
			}
			
		}

		return result / (4.*M_PI);
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
		
		double result = 0;

		for (unsigned i = 0; i < N; ++i)
		{
			auto const &c1 = corners.col(i);
			auto const &c2 = corners.col((i+1) % N);
			
			x_t dyxi = (c2-c1)/2.;
			
			for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
			{
				double xi = it->get_xi()(0);
				double w = it->get_w();
				
				x_t y = c1 * (1.-xi)/2. + c2 * (1.+xi)/2.;
				
				x_t rvec = y - x0;
				double r = rvec.norm();
				double z = -rvec(2);
				// square of lateral radius
				double R2 = rvec(0)*rvec(0) + rvec(1)*rvec(1);
				
				double integrand;
				/** \todo this decision could be made outside the loop over triangles */
				if (z == 0.0)
					integrand = 0.0;
				else
					integrand = z / std::abs(z) - z / r;
				
				double dtheta = (rvec(0) * dyxi(1) - rvec(1) * dyxi(0)) / R2;
				
				result += integrand * dtheta * w;
			} // loop over quadrature points
		} // end of loop over (triangles)

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
		
		double result = 0;

		for (unsigned i = 0; i < N; ++i)
		{
			auto const &c1 = corners.col(i);
			auto const &c2 = corners.col((i+1) % N);
			
			x_t dyxi = (c2-c1)/2.;
			
			for (auto it = quadrature_t::quadrature.begin(); it != quadrature_t::quadrature.end(); ++it)
			{
				double xi = it->get_xi()(0);
				double w = it->get_w();
				
				x_t y = c1 * (1.-xi)/2. + c2 * (1.+xi)/2.;
				
				x_t rvec = y - x0;
				double r = rvec.norm();
				double r2 = r * r;
				double z = -rvec(2);
				// square of lateral radius
				double R2 = rvec(0)*rvec(0) + rvec(1)*rvec(1);
				double R = std::sqrt(R2);
				
				double integrand;
				if (std::abs(z) > 1e-3)
				{
					double rdnx = -rvec.dot(nx) / r;
					integrand = -R2 / r2 / z * rdnx;
				}
				else
				{
					integrand = -nx(2) / R;
				}
    
				double dtheta = (rvec(0) * dyxi(1) - rvec(1) * dyxi(0)) / R2;
				
				result += integrand * dtheta * w;
			}
		}

		return result / (4.*M_PI);
	}
};




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
			test_field.get_elem().get_unit_normal());
		return result;
	}
};


} // end of namespace NiHu

#endif
