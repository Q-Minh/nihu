#include <iostream>

#include "library/guiggiani_1992.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"

unsigned const order = 5;

template <class elem_t>
double laplace_planar_analytic(elem_t const &elem, typename elem_t::x_t const &x0)
{
	typedef typename elem_t::x_t x_t;
	double res = 0.0;
	unsigned const N = elem_t::domain_t::num_corners;
	for (unsigned i = 0; i < N; ++i)
	{
		x_t c1 = elem.get_coords().col(i);
		x_t c2 = elem.get_coords().col((i + 1) % N);
		x_t d1 = c1 - x0, d2 = c2 - x0;
		x_t side = c2 - c1;
		x_t d = d1 - side.normalized() * side.normalized().dot(d1);

		double dd1 = d.normalized().dot(d1.normalized());
		double dd2 = d.normalized().dot(d2.normalized());

		res += (std::sqrt(1.0 - dd1*dd1) + std::sqrt(1.0 - dd2*dd2)) / d.norm();

//		res += (std::sin(phi2) + std::sin(phi1)) / d.norm();
	}
	return -res / (4.0 * M_PI);
}

void test_laplace_3d_quadratic(void)
{
	std::cout << "distortion test of laplace 3D kernel with quadratic triangle:\n";
	typedef tria_1_elem elem0_t;
	elem0_t::coords_t coords0;
	coords0 <<
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 0.0, 0.0;
	elem0_t elem0(coords0);

	typedef tria_2_elem elem_t;
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef trial_field_t test_field_t;
	typedef laplace_3d_HSP_kernel kernel_t;

	kernel_t kernel;

	elem_t::coords_t coords;
	coords <<
		0.0, 0.5, 1.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.5, 1.0, 0.5,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	trial_field_t::nset_t::shape_t I;
	for (double c = 1e-2; c < 1.0; c += 1e-2)
	{
		coords(0, 1) = coords(0, 3) = c;
		coords(1, 3) = coords(1,5) = 1.0 - c;
		elem_t elem(coords);
		guiggiani<test_field_t, trial_field_t, kernel_t, order> gui(elem, kernel);
		I.setZero();
		gui.integral(I);

		double I0 = laplace_3d_HSP_collocation_constant_triangle::eval(elem0, elem.get_center());

		std::cout << c << "\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;
	}
}


void test_laplace_3d_linear(void)
{
	std::cout << "distortion test of laplace 3D kernel with linear rectangle:\n";
	typedef quad_1_elem elem_t;
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t elem(coords);

	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef trial_field_t test_field_t;
	typedef laplace_3d_HSP_kernel kernel_t;

	kernel_t kernel;

	Eigen::Matrix<double, 1, 1> I;

	for (double c = 1e-2; c <= 2; c += 1e-2)
	{
		coords(1,2) = coords(1,3) = c;
		coords(0, 2) = 1.0 + c;
		coords(0, 3) = c;
		elem_t elem(coords);
		guiggiani<test_field_t, trial_field_t, kernel_t, order> gui(elem, kernel);
		I.setZero();
		gui.integral(I);

		double I0 = laplace_planar_analytic(elem, elem.get_center());

		std::cout << c << "\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << " " << I0 << " " << I << std::endl;
	}
}


void test_helmholtz_3d(void)
{
	typedef tria_1_elem elem0_t;
	typedef tria_2_elem elem_t;
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef trial_field_t test_field_t;
	typedef helmholtz_3d_HSP_kernel<double> kernel_t;

	elem0_t::coords_t coords0;
	coords0 <<
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 0.0, 0.0;
	coords0.row(1) *= 2.0;
	elem0_t elem0(coords0);

	elem_t::coords_t coords;
	coords <<
		0.0, 0.5, 1.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.5, 1.0, 0.5,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	coords.row(1) *= 2.0;
	elem_t elem(coords);

	kernel_t kernel(1.0);

	guiggiani<test_field_t, trial_field_t, kernel_t, order> gui(elem, kernel);

	Eigen::Matrix<std::complex<double>, 1, 1> I;
	I.setZero();
	gui.integral(I);
	std::complex<double> I0 = helmholtz_3d_HSP_collocation_constant_triangle::eval(
		elem0,
		elem.get_center(),
		1.0);

	std::cout << "I:\t" << I << std::endl;
	std::cout << "Ianal:\t" << I0 << std::endl;
	std::cout << "log10 error:\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;
}

int main(void)
{
//	test_laplace_3d_quadratic();
	test_laplace_3d_linear();
//	test_helmholtz_3d();

	return 0;
}
