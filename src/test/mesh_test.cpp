#include "mesh.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;

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

	std::for_each(
		mesh.begin<tria_1_elem>(),
		mesh.end<tria_1_elem>(),
		[] (tria_1_elem const &t) { std::cout << t.get_nodes() << std::endl; }
	);

	return 0;
}

