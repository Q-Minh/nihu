#include "rayleigh.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;
typedef tmp::unique<elem_vector>::type unique_vector;
typedef tria_1_elem::x_t x_t;

int main(void)
{
	double c[] = {
		0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 0
	};
	unsigned e[] = {
		4,  0, 1, 2, 3,
		3,  0, 1, 2, 0,
		3,  1, 2, 3, 0
	};

	Mesh<elem_vector> mesh;
	for (unsigned i = 0; i < 4; ++i)
		mesh.add_node(c+i*3);
	for (unsigned i = 0; i < 3; ++i)
		mesh.add_elem(e+i*5);

	field_points<x_t> field_p;
	field_p.add_point(x_t());

	rayleigh<elem_vector, constant_field> r(mesh, field_p, 1.0);
	r.eval();

	return 0;
}
