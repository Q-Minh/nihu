#include "core/inverse_mapping.hpp"
#include "library/lib_element.hpp"

#include <iostream>

int main(void)
{
	typedef NiHu::tria_2_elem elem_t;

	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 2.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0, 2.0, 1.0,
		0.0, 0.0, 1.0, 0.0, 0.0, 0.0;

	elem_t elem(coords);

	elem_t::x_t x;
	x << .2, .2, .3;

	NiHu::inverse_mapping<elem_t> im(elem);

	if (im.eval(x, 1e-5, 200))
	{
		std::cout << "Inverse mapping succeeded" << std::endl;
		std::cout << im.get_iter() << " iterations" << std::endl;
		std::cout << "Absolute error: " << im.get_error() << std::endl;
		std::cout << im.get_result() << std::endl;
	}
	else
	{
		std::cout << "Inverse mapping failed" << std::endl;
	}

	return 0;
}
