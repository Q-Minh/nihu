#include <boost/math/constants/constants.hpp>

#include <gtest/gtest.h>

#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/laplace_nearly_singular_integrals.hpp"
#include "nihu/library/laplace_singular_integrals.hpp"

#include "nihu/library/lib_element.hpp"
#include "nihu/core/field.hpp"
#include "nihu/library/line_1_gauss_shape_set.hpp"
#include "nihu/core/quadrature.hpp"


typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;

TEST(laplace_2d_elem_integrals, singular_2d_SLP_collocation)
{
	using namespace boost::math::double_constants;

	// check that analytical and general formulations yield the same result
	// for line_1 element with constant shape function
	typedef NiHu::line_1_elem elem_t;
	typedef NiHu::dirac_field<
		NiHu::field<elem_t, NiHu::line_1_gauss_shape_set>
	> gauss_field_t;
	
	elem_t::coords_t coords;
	coords <<
			-0.0, 1.0,
			0.0, 0.0;
	elem_t elem(coords);
	
	auto res1 = NiHu::laplace_2d_SLP_collocation_straight<
		gauss_field_t,
		NiHu::field<elem_t, NiHu::line_0_shape_set>
	>::eval(elem);
	
	auto res2 = NiHu::laplace_2d_SLP_collocation_curved<
		gauss_field_t,
		NiHu::field<elem_t, NiHu::line_0_shape_set>,
		1
	>::eval(elem);
	
	EXPECT_LE((res1-res2).norm() / res1.norm(), 1e-10);
	
	
	// check that analytical and general formulations yield the same result
	// for line_1 element with max. second order shape functions
	auto res3 = NiHu::laplace_2d_SLP_collocation_straight<
		gauss_field_t,
		NiHu::field<elem_t, NiHu::line_2_shape_set>
	>::eval(elem);
	
	auto res4 = NiHu::laplace_2d_SLP_collocation_curved<
		gauss_field_t,
		NiHu::field<elem_t, NiHu::line_2_shape_set>,
		1
	>::eval(elem);
	
	EXPECT_LE((res3-res4).norm() / res3.norm(), 1e-10);
	

	// check that analytical and general formulations yield the same result
	// for distorted quadratic but straight line element with constant shape function
	typedef NiHu::line_2_elem elem_2_t;
	typedef NiHu::dirac_field<
		NiHu::field<elem_2_t, NiHu::line_0_shape_set>
	> test_field_2_t;
	
	elem_2_t::coords_t coords2;
	coords2 <<
			0.0, 0.3, 1.0,
			0.0, 0.0, 0.0;
	elem_2_t elem2(coords2);
	
	double q = 1./two_pi * (.3 * (1-std::log(.3)) + .7 * (1.-std::log(.7)));
	
	auto res6 = NiHu::laplace_2d_SLP_collocation_curved<
		test_field_2_t,
		NiHu::field<elem_2_t, NiHu::line_0_shape_set>,
		10
	>::eval(elem2);
	
	EXPECT_LE(std::abs(q-res6(0,0)) / std::abs(res6(0,0)), 1e-10);
}


TEST(laplace_2d_elem_integrals, singular_2d_DLP_collocation)
{
	using namespace boost::math::double_constants;

	// check that analytical and general formulations yield the same result
	// for line_1 element with constant shape function
	typedef NiHu::line_1_elem elem_t;
	typedef NiHu::dirac_field<
		NiHu::field<elem_t, NiHu::line_0_shape_set>
	> test_field_t;

	elem_t::coords_t coords;
	coords <<
		0.0, 1.0,
		0.0, 0.0;
	elem_t elem(coords);

	auto res1 = NiHu::laplace_2d_DLP_collocation_curved<
		test_field_t,
		NiHu::field<elem_t, NiHu::line_0_shape_set>,
		10
	>::eval(elem);

	// check that analytical and general formulations yield the same result
	// for distorted quadratic but straight line element with constant shape function
	typedef NiHu::line_2_elem elem_2_t;

	elem_2_t::coords_t coords2;
	coords2 <<
		0.0, 0.3, 1.0,
		0.0, 0.0, 0.0;
	elem_2_t elem2(coords2);

	auto res2 = NiHu::laplace_2d_DLP_collocation_curved<
		test_field_t,
		NiHu::field<elem_2_t, NiHu::line_0_shape_set>,
		10
	>::eval(elem2);

 //	EXPECT_LE((res1 - res2).norm() / res2.norm(), 1e-10);
}


TEST(laplace_2d_elem_integrals, singular_2d_HSP_collocation)
{
	double I1, I2;
	
	typedef NiHu::line_1_gauss_shape_set test_shape_t;
	
	{
		typedef NiHu::line_1_elem elem_t;
		typedef NiHu::line_1_gauss_shape_set trial_shape_set_t;
		
		typedef NiHu::field<elem_t, test_shape_t> base_field_t;
		typedef NiHu::dirac_field<base_field_t> test_field_t;
		typedef NiHu::field<elem_t, trial_shape_set_t> trial_field_t;

		// create element
		elem_t::coords_t coords;
		coords <<
			-0.0, 1.0,
			0.0, 0.0;
		elem_t elem(coords);
		
		auto res = NiHu::laplace_2d_HSP_collocation_straight<
			test_field_t, trial_field_t
			>::eval(elem);
		I1 = res.row(0).sum();
	}
	
	{
		typedef NiHu::line_2_elem elem_t;
		typedef NiHu::line_2_shape_set trial_shape_t;
		
		typedef NiHu::field<elem_t, test_shape_t> base_field_t;
		typedef NiHu::dirac_field<base_field_t> test_field_t;
		typedef NiHu::field<elem_t, trial_shape_t> trial_field_t;

		// create element
		elem_t::coords_t coords;
		coords <<
			0.0, 0.5, 1.0,
			0.0, 0.0, 0.0;
		elem_t elem(coords);
		
		auto res = NiHu::laplace_2d_HSP_collocation_curved<
			test_field_t, trial_field_t, 20
			>::eval(elem);
		I2 = res.row(0).sum();
	}
	
	EXPECT_LE(std::abs(I1-I2)/std::abs(I1), 1e-10);
}


TEST(laplace_2d_elem_integrals, singular_2d_SLP_galerkin)
{
	using namespace boost::math::double_constants;

	// check that analytical and general formulations yield the same result
	// for line_1 element with constant shape function
	typedef NiHu::line_1_elem elem_t;
	typedef NiHu::field<elem_t, NiHu::line_1_shape_set> test_field_t;
	typedef NiHu::field<elem_t, NiHu::line_1_shape_set> trial_field_t;

	elem_t::coords_t coords;
	coords <<
		0.0, 1.0,
		0.0, 0.0;
	elem_t elem(coords);

	auto res1 = NiHu::laplace_2d_SLP_galerkin_face_general<
		test_field_t, trial_field_t, 10
	>::eval(elem);

	double d1, d2;
	NiHu::laplace_2d_SLP_galerkin_face_linear_line::eval(elem, d1, d2);

	std::cout << res1 << std::endl;
	std::cout << d1 << ' ' << d2 << std::endl;
}
