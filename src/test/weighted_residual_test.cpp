#include "../bem/weighted_residual.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;
typedef Mesh<elem_vector> mesh_t;
typedef function_space<mesh_t, isoparametric_field> func_space_t;

#include <iostream>

int main(void)
{
	double c[] = {
		0, 0, 0,
		1, 0, 0,
		2, 0, 0,
		0, 1, 0,
		1, 1, 0,
		2, 1, 0,
	};
	unsigned e[] = {
		4,  0, 1, 4, 3,
		4,  1, 2, 5, 4,
		3,  1, 2, 5, 0
	};

	Mesh<elem_vector> mesh;
	for (unsigned i = 0; i < 6; ++i)
		mesh.add_node(c+i*3);
	for (unsigned i = 0; i < 3; ++i)
		mesh.add_elem(e+i*5);

	func_space_t func(mesh);
	weighted_residual<green_G_kernel, func_space_t, func_space_t>wr(func, func);

	int a;
	wr.eval(a);

	return 0;
}

