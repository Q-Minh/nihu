#include "util/brick.hpp"
#include "library/location.hpp"

int main(void)
{
	typedef build<
		location<space_3d>::brick,
		build<area_elem<space_3d>::brick>
	> input_t;

	typedef merge<input_t, build<area_elem<space_3d>::brick> >::type input2;

	typename quad_1_elem::coords_t c;
	c << 0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;
	quad_1_elem q(c);

	input2 input(q, quad_1_elem::domain_t::get_center());
	std::cout << input.get_x() << std::endl << input.get_dA() << std::endl;
}

