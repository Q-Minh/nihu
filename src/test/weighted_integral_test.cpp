#include "../bem/weighted_integral.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;
typedef Mesh<elem_vector> mesh_t;

#include <iostream>

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

	function_space<mesh_t, isoparametric_field> func(mesh);
	weighted_integral<function_space<Mesh<elem_vector>, isoparametric_field>, green_G_kernel> wi(func);

	wi.eval();

	return 0;
}

