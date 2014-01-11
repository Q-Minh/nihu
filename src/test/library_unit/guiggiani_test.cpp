#include <iostream>
#include "library/guiggiani_1992.hpp"

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

void test_3d(void)
{
	typedef quad_1_elem elem_t;
	typedef elem_t::xi_t xi_t;
	typedef field_view<elem_t, field_option::constant> field_t;
	typedef laplace_3d_HSP_kernel kernel_t;

	elem_t::coords_t coords;
	coords <<
		-1.5, +1.5, +1.5, -1.5,
		-1.5, -1.5, +1.5, +1.5,
		0, 0, 0, 0;
	/*
	coords <<
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 0.0;
	*/
	elem_t elem(coords);
	xi_t xi0 = elem_t::domain_t::get_center();
	kernel_t kernel;

	guiggiani<field_t, kernel_t> gui(elem, kernel);

	field_t::nset_t::shape_t I;
	I.setZero();
	gui.integral(xi0, I);
	double I0 = anal_3d(elem, xi0);

	std::cout << "I:\t" << I << std::endl;
	std::cout << "Ianal:\t" << I0 << std::endl;
	std::cout << "log10 error:\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;
}


void test_2d(void)
{
	typedef line_2_elem elem_t;
	typedef elem_t::xi_t xi_t;
	typedef field_view<elem_t, field_option::constant> field_t;
	typedef laplace_2d_HSP_kernel kernel_t;

	elem_t::coords_t coords;
	coords <<
		-1.0, 0.5, +1.0,
		0, 0, 0;
	elem_t elem(coords);
	xi_t xi0 = elem_t::domain_t::get_center();
	kernel_t kernel;

	guiggiani<field_t, kernel_t> gui(elem, kernel);

	field_t::nset_t::shape_t I;
	I.setZero();
	gui.integral(xi0, I);
//	double I0 = anal_2d(elem, xi0);

	std::cout << "I:\t" << I << std::endl;
//	std::cout << "Ianal:\t" << I0 << std::endl;
//	std::cout << "log10 error:\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;
}


int main(void)
{
	test_2d();

	return 0;
}
