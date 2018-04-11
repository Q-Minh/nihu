#include <gtest/gtest.h>

#include "library/lib_element.hpp"
#include "core/field.hpp"
#include "library/helmholtz_kernel.hpp"
#include "core/double_integral.hpp"
#include "../benchmark/EAA/radiatterer/bem/quad_1_gauss_elem.hpp"

#include "library/helmholtz_singular_integrals.hpp"


typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cVector;

TEST(HelmholtzElemIntegrals, Regular)
{
	typedef NiHu::quad_1_elem elem_t;

	// create test_elem
	elem_t::coords_t test_coords;
	test_coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t test_elem(test_coords);
	std::cout << test_elem.get_id() << std::endl;

	// create trial_elem
	elem_t::coords_t trial_coords;
	trial_coords <<
		5.0, 6.0, 6.0, 5.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t trial_elem(trial_coords);
	std::cout << trial_elem.get_id() << std::endl;

	// create test_field
	typedef NiHu::field_view<elem_t, NiHu::field_option::constant> base_field_t;
	typedef NiHu::dirac_field<base_field_t> test_field_t;
	test_field_t const &test_field = NiHu::dirac(NiHu::constant_view(test_elem));

	// instantiate kernel
	typedef NiHu::helmholtz_3d_HSP_kernel<double> kernel_t;
	double k = 2;
	kernel_t kernel(k);

	{
		// create trial_field
		typedef NiHu::quad_2_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, false>());
		auto I = res * cVector::Constant(9, 1, 1.0);

		std::cout << I << std::endl;
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
		auto I = res * cVector::Constant(4, 1, 1.0);

		std::cout << I << std::endl;
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
		auto I = res * cVector::Constant(1, 1, 1.0);

		std::cout << I << std::endl;
	}

}



TEST(HelmholtzElemIntegrals, Singular)
{
	typedef NiHu::quad_1_elem elem_t;

	// create two overlapping elements
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 1.0, 0.1,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t test_elem(coords), trial_elem(coords);

	// create test_field
	typedef NiHu::field_view<elem_t, NiHu::field_option::constant> base_field_t;
	typedef NiHu::dirac_field<base_field_t> test_field_t;
	test_field_t const &test_field = NiHu::dirac(NiHu::constant_view(test_elem));

	// instantiate kernel
	typedef NiHu::helmholtz_3d_HSP_kernel<double> kernel_t;
	double k = 2;
	kernel_t kernel(k);

	auto x0 = test_elem.get_center();
	std::complex<double> I0 = NiHu::helmholtz_3d_HSP_collocation_constant_plane<15>::eval(test_elem, x0, k);

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
		std::complex<double> I = (res * cVector::Constant(9, 1, 1.0))(0);
		double err = std::abs(I-I0) / std::abs(I0);
		EXPECT_LE(err, eps);

	}

	{
		// create trial_field
		typedef quad_1_gauss_shape_set trial_shape_set;
		typedef NiHu::field<elem_t, trial_shape_set> trial_field_t;
		trial_field_t::dofs_t dofs;
		trial_field_t trial_field(trial_elem, dofs);

		auto res =  NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
			kernel, test_field, trial_field,
			std::integral_constant<bool, true>());
		std::complex<double> I = (res * cVector::Constant(4, 1, 1.0))(0);
		double err = std::abs(I-I0) / std::abs(I0);
		EXPECT_LE(err, eps);
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
		std::complex<double> I = (res * cVector::Constant(4, 1, 1.0))(0);
		double err = std::abs(I-I0) / std::abs(I0);
		EXPECT_LE(err, eps);
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
		std::complex<double> I = (res * cVector::Constant(1, 1, 1.0))(0);
		double err = std::abs(I-I0) / std::abs(I0);
		EXPECT_LE(err, eps);
	}
}

