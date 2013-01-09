#include <iostream>
#include "element.hpp"

int main(void)
{
	Matrix<double,4,3> coords;
	coords <<
		0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 0;
	if (shape_set_converter<quad_1_shape_set, parallelogram_shape_set>::eval(coords))
	{
		parallelogram_elem e(coords.topRows(3));
		std::cout << coords << ": parallelogram" << std::endl;
	}
	else
	{
		quad_1_elem e(coords);
		std::cout << coords << ": general quad_1_element" << std::endl;
	}

	return 0;
}

