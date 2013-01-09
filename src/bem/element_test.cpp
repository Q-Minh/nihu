#include <iostream>
#include "element.hpp"

int main(void)
{
	Matrix<double,4,3> coords2;
	coords2 <<
		0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 1;
	if (lset_converter<quad_1_lset, tria_1_lset>::eval(coords2))
		parallelogram_elem e(coords2.topRows(3));

	return 0;
}

