#include "util/brick.hpp"
#include "library/location.hpp"

int main(void)
{
	typedef
	build<
		normal_jacobian<space_3d>::brick,
		build<
			location<space_3d>::brick
		>
	> input_t;
	
	typedef merge<
		input_t,
		build<normal_jacobian<space_3d>::brick>
	>::type input_opt;

	typename quad_1_elem::coords_t c;
	c << 0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;
	quad_1_elem q(c);

	input_opt input(q, quad_1_elem::domain_t::get_center());
	std::cout << input.get_x() << std::endl
		<< input.get_normal() << std::endl
		<< input.get_jacobian() << std::endl;
}

