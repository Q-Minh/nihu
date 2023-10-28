#include <gtest/gtest.h>

#include <boost/math/constants/constants.hpp>

#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/laplace_nearly_singular_integrals.hpp"
#include "nihu/library/laplace_singular_integrals.hpp"

#include "nihu/library/lib_element.hpp"
#include "nihu/core/field.hpp"
#include "nihu/core/double_integral.hpp"
#include "nihu/library/quad_1_gauss_shape_set.hpp"


typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;

TEST(LaplaceElemIntegrals, Regular)
{
	typedef NiHu::tria_1_elem elem_t;
	
	elem_t::coords_t coords;
	coords << 0.0, 1.0, 0.0,
			0.0, 0.0, 0.7,
			0.0, 0.0, 0.0;
	elem_t elem(coords);
			
	elem_t::x_t x0;
	x0 << .3, .3, .3;
	
	double g = NiHu::laplace_3d_SLP_collocation_constant_plane_nearly_singular::eval(elem, x0);
	
	// from Matlab with 160-th order quadrature
	double g0 = 0.070405172409714;
	
	double rel_error = std::abs(g0 - g) / std::abs(g0);
	EXPECT_LE(rel_error, 1e-6);
}


static double helper(double x, double t0, double a, double b, double const *c)
{
	return (c[0] * cos(a + x)) / b 
		+ 2. * ( c[1] * sin(a + t0) + c[2] * cos(a + t0)) * atanh(cos(a) - sin(a) * tan(x/2.)) 
		- (c[1] * sin(t0 - x) + c[2] * cos(t0 - x)) * (log(b / sin(a + x)) + 1.) 
		+ c[3] / 2. * b * (sin(2. * (a + t0)) * log(1./tan((a + x)/2.)) - 2. * sin(a + 2. * t0 - x));
}


TEST(LaplaceElemIntegrals, Singular_HSP_Linear_NSet)
{
	using namespace boost::math::double_constants;

	typedef NiHu::quad_1_elem elem_t;
	
	typedef NiHu::quad_1_gauss_shape_set test_shape_set_t;
	typedef NiHu::field<elem_t, test_shape_set_t> base_field_t;
	typedef NiHu::dirac_field<base_field_t> test_field_t;
	
	typedef NiHu::quad_1_gauss_shape_set trial_shape_set_t;
	typedef NiHu::field<elem_t, trial_shape_set_t> trial_field_t;
	
	typedef NiHu::laplace_3d_HSP_kernel kernel_t;
	
	typedef NiHu::double_integral<kernel_t, test_field_t, trial_field_t> integrator_t;

	// create two overlapping elements
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t elem(coords);
	
	// create test field
	base_field_t base_field(elem, base_field_t::dofs_t());
	test_field_t const &test_field = NiHu::dirac(base_field);
	
	// create trial_field
	trial_field_t trial_field(elem, trial_field_t::dofs_t());

	// instantiate kernel
	kernel_t kernel;

	auto res = integrator_t::eval(
		kernel, test_field, trial_field, std::true_type());
	std::cout << res.row(0) << std::endl;
	
	
	// compute analytically
	
	// location of collocation point in physical coordinates
	double q = .5 / sqrt(3.);
	Eigen::Matrix<double, 3, 1> x0;
	x0 << .5-q, .5-q, 0.;
	
	// divide the element into triangles
	double R[4], theta[4], alpha[4], t0[4];
	plane_element_helper(elem, x0, R, theta, alpha);
	for (int i = 0; i < 4; ++i)
		t0[i] = atan2(coords(1,i)-x0(0), coords(0,i)-x0(0));
	
	// shape function coefficients
	double coeffs[4][4] = {
	// const,      x,        y,              xy
		{1, -1/(2*q), -1/(2*q),  1/(2*q)/(2*q)},
		{0,  1/(2*q),        0, -1/(2*q)/(2*q)},
		{0,        0,  1/(2*q), -1/(2*q)/(2*q)},
		{0,        0,        0,  1/(2*q)/(2*q)}
	};
	
	Eigen::Matrix<double, 1, 4> res0;
	res0.setZero();
		
	for (int tr = 0; tr < 4; ++tr)
	{
		double b = R[tr] * sin(alpha[tr]);
		for (int c = 0; c < 4; ++c)
			res0(c) += helper(theta[tr], -t0[tr], alpha[tr], b, coeffs[c])
				- helper(0, -t0[tr], alpha[tr], b, coeffs[c]);
	}
	res0 /= (4. * pi);
	
	// compare numerical and analytical solutions
	double rel_error = (res0 - res.row(0)).norm() / res0.norm();
	EXPECT_LE(rel_error, 1e-4);
}
