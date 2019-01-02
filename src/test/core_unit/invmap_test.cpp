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

	auto res = NiHu::inverse_mapping<elem_t>::eval(elem, x, 200, 1e-5);

	elem_t::x_t y = elem.get_x(res.topRows(2));
	elem_t::x_t n = elem.get_normal(res.topRows(2));
	elem_t::x_t x0 = y + n * res(2);

	return 0;
}
