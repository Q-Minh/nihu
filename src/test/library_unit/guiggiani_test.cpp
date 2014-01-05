#include <iostream>
#include "library/guiggiani_1992.hpp"

typedef tria_1_elem elem_t;
typedef elem_t::xi_t xi_t;
typedef field_view<elem_t, field_option::constant> field_t;
typedef laplace_3d_HSP_kernel kernel_t;

double anal(elem_t const &elem, xi_t const &xi0)
{
	typedef elem_t::x_t x_t;
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

int main(void)
{
	elem_t::coords_t coords;
	/*
	coords <<
		-1.5, +1.5, +1.5, 0.0,
		-1.5, -1.5, +1.5, +1.5,
		0, 0, 0, 0;
		*/
	coords <<
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 0.0, 0.0;
	elem_t elem(coords);
	elem_t::xi_t xi0 = elem_t::domain_t::get_center();

	guiggiani<field_t, kernel_t> gui(elem);

	field_t::nset_t::shape_t I;
	I.setZero();
	gui.integral(xi0, I);
	double I0 = anal(elem, xi0);

	std::cout << "I:\t" << I << std::endl;
	std::cout << "Ianal:\t" << I0 << std::endl;
	std::cout << "log10 error:\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;

	return 0;
}
