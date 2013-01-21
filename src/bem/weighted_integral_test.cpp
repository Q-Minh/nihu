#include "weighted_integral.hpp"

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
    for (unsigned n = 0; n < 1000; ++n)
    {
        for (unsigned i = 0; i < 4; ++i)
            mesh.add_node(c+i*3);
        for (unsigned i = 0; i < 3; ++i)
            mesh.add_elem(e+i*5);
    }

	function_space<Mesh<elem_vector>, constant_field> con_func(mesh);
	weighted_integral<function_space<Mesh<elem_vector>, constant_field>, green_kernel> wi(con_func);

	gauss_quad<tria_domain, 5>::init();

	wi.do_it();

	return 0;
}

