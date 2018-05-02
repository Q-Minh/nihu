#include <gtest/gtest.h>

#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "library/laplace_hfp.hpp"

#include "library/lib_element.hpp"
#include "core/field.hpp"
#include "library/line_1_gauss_shape_set.hpp"
#include "core/quadrature.hpp"


typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;

TEST(LaplaceElemIntegrals, Singular_2d_SLP)
{
	typedef NiHu::line_1_elem elem_t;
	typedef NiHu::line_1_gauss_shape_set test_shape_t;
	typedef NiHu::field<elem_t, test_shape_t> base_field_t;
	typedef NiHu::dirac_field<base_field_t> test_field_t;
	
	typedef NiHu::line_1_gauss_shape_set trial_shape_set_t;
	typedef NiHu::field<elem_t, trial_shape_set_t> trial_field_t;

	// create two overlapping elements
	elem_t::coords_t coords;
	coords <<
		-0.0, 1.0,
		0.0, 0.0;
	elem_t elem(coords);
	
	auto res = NiHu::laplace_2d_SLP_collocation_line<test_field_t, trial_field_t>::eval(elem);
	std::cout << res << std::endl;

	
	res = NiHu::laplace_2d_HSP_collocation_line<test_field_t, trial_field_t>::eval(elem);
	std::cout << res << std::endl;

	double d = NiHu::laplace_2d_HSP_collocation_constant_line::eval(elem, elem.get_x(test_shape_t::corner_at(0)));
	std::cout << d << std::endl;
}

TEST(LaplaceElemIntegrals, Singular_2d_HSP)
{
	typedef NiHu::line_2_elem elem_t;
	typedef NiHu::line_1_gauss_shape_set test_shape_t;
	typedef NiHu::field<elem_t, test_shape_t> base_field_t;
	typedef NiHu::dirac_field<base_field_t> test_field_t;
	typedef NiHu::line_2_shape_set trial_shape_t;
	typedef NiHu::field<elem_t, trial_shape_t> trial_field_t;

	// create element
	elem_t::coords_t coords;
	coords <<
		-0.0, 0.5, 1.0,
		0.0, 0.0, 0.0;
	elem_t elem(coords);
	
	Eigen::Matrix<double, 2, 3> result_slp;
	result_slp.setZero();
	NiHu::laplace_2d_SLP_collocation_general<test_field_t, trial_field_t>::eval(elem, result_slp);
	
	std::cout << "SLP: " << std::endl << result_slp << std::endl;
	
	Eigen::Matrix<double, 2, 3> result_hsp;
	result_hsp.setZero();
	NiHu::laplace_2d_HSP_collocation_general<test_field_t, trial_field_t>::eval(elem, result_hsp);
	
	std::cout << "HSP: " << std::endl << result_hsp << std::endl;
}

