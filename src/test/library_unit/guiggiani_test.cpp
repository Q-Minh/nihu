#include <iostream>
#include "library/guiggiani_1992.hpp"

typedef quad_1_elem elem_t;
typedef field_view<quad_1_elem, field_option::constant> field_t;
typedef laplace_3d_HSP_kernel kernel_t;


int main(void)
{
	elem_t::coords_t coords;
	coords <<
		-1, +1, +1, -1,
		-1, -1, +1, +1,
		1, 0, 0, 0;
	elem_t elem(coords);

	typedef guiggiani<field_t, kernel_t> gui_t;

	gui_t gui(elem);

	auto xi0 = elem_t::domain_t::get_center();

	field_t::nset_t::shape_t I1, I2;
	I1.setZero();
	I2.setZero();
	gui.line_integrals(xi0, I1, I2);

	std::cout << "I1: " << I1 << std::endl;
	std::cout << "I2: " << I2 << std::endl;

	/*
	std::cout << "4 pi Fm2: " << guiggiani_t::Fm2(elem, xi0, theta) * 4 * M_PI << std::endl;
	std::cout << "4 pi Fm1: " << guiggiani_t::Fm1(elem, xi0, theta) * 4 * M_PI << std::endl;
	*/
	return 0;
}
