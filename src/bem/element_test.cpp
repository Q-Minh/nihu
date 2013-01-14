#include <iostream>
#include "element.hpp"

int main(void)
{
	Matrix<double,4,3> coords;
	coords <<
		0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 1;
	if (shape_set_converter<quad_1_shape_set, parallelogram_shape_set>::eval(coords))
	{
		parallelogram_elem::coords_type c;
		c.row(0) = coords.row(0);
		c.row(1) = coords.row(1);
		c.row(2) = coords.row(3);
		parallelogram_elem e(c);

		std::cout << coords << ": parallelogram" << std::endl;
		std::cout << e.get_dx(Matrix<double, 2, 1>::Zero()) << std::endl;
	}
	else
	{
		quad_1_elem e(coords);
		std::cout << coords << ": general quad_1_element" << std::endl;
		std::cout << e.get_dx(Matrix<double, 2, 1>::Zero()) << std::endl;
	}

	return 0;
}

