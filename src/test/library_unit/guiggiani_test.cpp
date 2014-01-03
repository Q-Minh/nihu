#include <iostream>
#include "library/guiggiani_1992.hpp"

typedef quad_1_elem elem_t;
typedef elem_t::xi_t xi_t;
typedef field_view<quad_1_elem, field_option::constant> field_t;
typedef laplace_3d_HSP_kernel kernel_t;

double anal(elem_t const &elem, xi_t const &xi0)
{
	double res = 0;
	elem_t::x_t x0 = elem.get_x(xi0);
	unsigned const N = elem_t::domain_t::num_corners;
	for (unsigned i = 0; i < N; ++i)
	{
		elem_t::x_t c1 = elem.get_coords().col(i);
		elem_t::x_t c2 = elem.get_coords().col((i + 1) % N);
		elem_t::x_t d1 = c1 - x0;
		elem_t::x_t d2 = c2 - x0;
		elem_t::x_t side = c2 - c1;
		elem_t::x_t d = d1 - side.normalized() * side.normalized().dot(d1);

		double phi1 = std::acos(d.normalized().dot(d1.normalized()));
		double phi2 = std::acos(d.normalized().dot(d2.normalized()));

		res += d.norm() * (std::sin(phi2) + std::sin(phi1));
	}
	return -res / (4.0 * M_PI);
}

int main(void)
{
	elem_t::coords_t coords;
	coords <<
		-.5, +.5, +.5, -.5,
		-.5, -.5, +.5, +.5,
		0, 0, 0, 0;
	elem_t elem(coords);

	typedef guiggiani<field_t, kernel_t> gui_t;

	gui_t gui(elem);

	auto xi0 = elem_t::domain_t::get_center();

	field_t::nset_t::shape_t I;
	I.setZero();
	gui.integral(xi0, I);

	std::cout << "I: " << I << std::endl;

	std::cout << "Ianal:" << anal(elem, xi0) << std::endl;

	/*
	std::cout << "4 pi Fm2: " << guiggiani_t::Fm2(elem, xi0, theta) * 4 * M_PI << std::endl;
	std::cout << "4 pi Fm1: " << guiggiani_t::Fm1(elem, xi0, theta) * 4 * M_PI << std::endl;
	*/
	return 0;
}
