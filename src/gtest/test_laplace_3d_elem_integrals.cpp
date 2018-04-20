#include <gtest/gtest.h>

#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

#include "library/lib_element.hpp"
#include "core/field.hpp"
#include "core/double_integral.hpp"
#include "library/quad_1_gauss_field.hpp"


typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;

TEST(LaplaceElemIntegrals, Regular)
{
	typedef NiHu::quad_1_elem elem_t;

	// create test_elem
	elem_t::coords_t test_coords;
	test_coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t test_elem(test_coords);

	// create trial_elem
	elem_t::coords_t trial_coords;
	trial_coords <<
		5.0, 6.0, 6.0, 5.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t trial_elem(trial_coords);

	// create test_field
	typedef NiHu::field_view<elem_t, NiHu::field_option::constant> base_field_t;
	typedef NiHu::dirac_field<base_field_t> test_field_t;
	test_field_t const &test_field = NiHu::dirac(NiHu::constant_view(test_elem));

	// instantiate kernel
	typedef NiHu::laplace_3d_HSP_kernel kernel_t;
	kernel_t kernel;

	{
		// create trial_field
		typedef NiHu::quad_2_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, false>());
		auto I = res * dVector::Constant(9, 1, 1.0);
	}

	{
		// create trial_field
		typedef NiHu::quad_1_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, false>());
		auto I = res * dVector::Constant(4, 1, 1.0);
	}


	{
		// create trial_field
		typedef NiHu::quad_0_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, false>());
		auto I = res * dVector::Constant(1, 1, 1.0);
	}

}



TEST(LaplaceElemIntegrals, Singular)
{
	typedef NiHu::quad_1_elem elem_t;

	// create two overlapping elements
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t test_elem(coords), trial_elem(coords);

	// create test_field
	typedef NiHu::quad_1_gauss_shape_set test_shape_set;
	typedef NiHu::field<elem_t, test_shape_set> base_field_t;
	typedef NiHu::dirac_field<base_field_t> test_field_t;
		base_field_t::dofs_t dofs;
		base_field_t base_field(test_elem, dofs);
	
	test_field_t const &test_field = NiHu::dirac(base_field);

	// instantiate kernel
	typedef NiHu::laplace_3d_HSP_kernel kernel_t;
	kernel_t kernel;

	auto x0 = test_elem.get_center();
	double I0 = NiHu::laplace_3d_HSP_collocation_constant_plane::eval(test_elem, x0);

	double eps = 1e-3;


	{
		// create trial_field
		typedef NiHu::quad_2_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, true>());
		//double I = (res * dVector::Constant(9, 1, 1.0))(0);
		//double err = std::abs(I-I0) / std::abs(I0);
		//EXPECT_LE(err, eps);

	}

	{
		// create trial_field
		typedef NiHu::quad_1_gauss_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, true>());
		std::cout << res << std::endl;
		std::cout << (res * dVector::Constant(4, 1, 1.0)) << std::endl;
		//double err = std::abs(I-I0) / std::abs(I0);
		//EXPECT_LE(err, eps);
	}

	{
		// create trial_field
		typedef NiHu::quad_1_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, true>());
		//double I = (res * dVector::Constant(4, 1, 1.0))(0);
		//double err = std::abs(I-I0) / std::abs(I0);
		//EXPECT_LE(err, eps);
	}


	{
		// create trial_field
		typedef NiHu::quad_0_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, true>());
		//double I = (res * dVector::Constant(1, 1, 1.0))(0);
		//double err = std::abs(I-I0) / std::abs(I0);
		//EXPECT_LE(err, eps);
	}
}

