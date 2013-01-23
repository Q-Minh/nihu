#include "rayleigh.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef rayleigh<elem_type_vector, isoparametric_field> rayleigh_t;


int main(void)
{
	double c[] = {
		0., 0., 0.,
		1., 0., 0.,
		2., 0., 0.,
		0., 1., 0.,
		1., 1., 0.,
		2., 1., 0.,
		0., 2., 0.,
		1., 2., 0.,
		2., 2., 0.
	};
	unsigned e[] = {
		4,  0, 1, 4, 3,
		4,  1, 2, 5, 4,
		4,  3, 4, 7, 6,
		3,  4, 5, 8, 0,
		3,  4, 8, 7, 0
	};

	Mesh<elem_type_vector> mesh;
	for (unsigned i = 0; i < sizeof(c)/sizeof(c[0])/3; ++i)
		mesh.add_node(c+i*3);
	for (unsigned i = 0; i < sizeof(e)/sizeof(e[0])/5; ++i)
		mesh.add_elem(e+i*5);

	rayleigh_t r(mesh);
	auto result = r.eval(rayleigh_t::x_t::Ones(), 1.0);
	std::cout << "result : " << std::endl << result << std::endl;

	return 0;
}

