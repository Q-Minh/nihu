#include <iostream>
#include "element.hpp"
#include "elem_descriptor.hpp"

int main(void)
{
	/*
	Matrix<double,4,3> coords;
	coords <<
		0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 1;
	if (shape_set_converter<quad_1_shape_set, parallelogram_shape_set>::eval<3>(coords))
	{
		std::cout << coords << ": parallelogram" << std::endl;

		parallelogram_elem::coords_type c;
		c.row(0) = coords.row(0);
		c.row(1) = coords.row(1);
		c.row(2) = coords.row(3);
		parallelogram_elem e(c);

		parallelogram_elem::xi_type xi = parallelogram_elem::xi_type::Zero();

		std::cout << "normal: " << e.get_normal(xi) << std::endl;
		std::cout << e.get_dx(xi) << std::endl;
	}
	else
	{
		std::cout << coords << ": general quad_1_element" << std::endl;

		quad_1_elem e(coords);

		quad_1_elem::xi_type xi = parallelogram_elem::xi_type::Zero();

		std::cout << "normal: " << e.get_normal(xi) << std::endl;
		std::cout << e.get_dx(xi) << std::endl;
	}
	*/

	gauss_quad<line_domain, 5>::init();

	Matrix<double,2,2> line_coords;
	line_coords <<
		0, 0,
		1, 1;
	line_1_elem l(line_coords);
	ElemAccelerator<line_1_elem, 5> a(l);

	for (auto i = a.begin(); i != a.end(); ++i)
		std::cout << i->get_x() << std::endl;

	return 0;
}

