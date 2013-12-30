#include <iostream>
#include "library/guiggiani_1992.hpp"

typedef quad_1_elem elem_t;
typedef field_view<quad_1_elem, field_option::isoparametric> field_t;

typedef guiggiani_hypersingular<laplace_3d_HSP_kernel, field_t> guiggiani_t;

int main(void)
{
	elem_t::coords_t coords;
	coords <<
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;
	elem_t elem(coords);
	auto xi0 = elem_t::domain_t::get_corner(1);

	double theta = 0.0;
	std::cout << guiggiani_t::Fm2(elem, xi0, theta) * 4 * M_PI << std::endl;
	std::cout << guiggiani_t::Fm1(elem, xi0, theta) * 4 * M_PI << std::endl;
	return 0;
}
