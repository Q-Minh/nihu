#include <iostream>

#include "library/guiggiani_1992.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"

template <class elem_t>
double anal_3d(elem_t const &elem, typename elem_t::xi_t const &xi0)
{
	typedef typename elem_t::x_t x_t;
	double res = 0;
	x_t x0 = elem.get_x(xi0);
	unsigned const N = elem_t::domain_t::num_corners;
	for (unsigned i = 0; i < N; ++i)
	{
		x_t c1 = elem.get_coords().col(i);
		x_t c2 = elem.get_coords().col((i + 1) % N);
		x_t d1 = c1 - x0, d2 = c2 - x0;
		x_t side = c2 - c1;
		x_t d = d1 - side.normalized() * side.normalized().dot(d1);

		double phi1 = std::acos(d.normalized().dot(d1.normalized()));
		double phi2 = std::acos(d.normalized().dot(d2.normalized()));

		res += (std::sin(phi2) + std::sin(phi1)) / d.norm();
	}
	return -res / (4.0 * M_PI);
}

void test_laplace_3d(void)
{
	typedef tria_1_elem elem0_t;
	typedef tria_2_elem elem_t;
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef trial_field_t test_field_t;
	typedef laplace_3d_HSP_kernel kernel_t;

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

	kernel_t kernel;

	guiggiani<test_field_t, trial_field_t, kernel_t> gui(elem, kernel);

	trial_field_t::nset_t::shape_t I;
	I.setZero();
	gui.integral(I);

	double I0 = laplace_3d_HSP_collocation_constant_triangle::eval(elem0);

	std::cout << "I:\t" << I << std::endl;
	std::cout << "Ianal:\t" << I0 << std::endl;
	std::cout << "log10 error:\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;
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

	guiggiani<test_field_t, trial_field_t, kernel_t> gui(elem, kernel);

	Eigen::Matrix<std::complex<double>, 1, 1> I;
	I.setZero();
	gui.integral(I);
	std::complex<double> I0 = helmholtz_3d_HSP_collocation_constant_triangle::eval(elem0, 1.0);

	std::cout << "I:\t" << I << std::endl;
	std::cout << "Ianal:\t" << I0 << std::endl;
	std::cout << "log10 error:\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;
}


int main(void)
{
	test_laplace_3d();
	test_helmholtz_3d();

	return 0;
}
