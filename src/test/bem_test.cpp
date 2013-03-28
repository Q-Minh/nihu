#include "../bem/bem.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef Mesh<elem_type_vector> mesh_t;
typedef mesh_t::x_t x_t;
typedef bem_radiation<green_G_kernel, elem_type_vector, isoparametric_field> iso_bem_t;
typedef bem_radiation<green_HG_kernel, elem_type_vector, constant_field> const_bem_t;

int main(void)
{
	double c[][3] = {
		{0., 0., 0.},
		{1., 0., 0.},
		{2., 0., 0.},
		{0., 1., 0.},
		{1., 1., 0.},
		{2., 1., 0.},
		{0., 2., 0.},
		{1., 2., 0.},
		{2., 2., 0.}
	};
	unsigned e[][5] = {
		{4,  0, 1, 4, 3},
		{4,  1, 2, 5, 4},
		{4,  3, 4, 7, 6},
		{3,  4, 5, 8, 0},
		{3,  4, 8, 7, 0}
	};

	mesh_t mesh;
	for (unsigned i = 0; i < sizeof(c)/sizeof(c[0]); ++i)
		mesh.add_node(c[i]);
	for (unsigned i = 0; i < sizeof(e)/sizeof(e[0]); ++i)
		mesh.add_elem(e[i]);

	x_t x0 = x_t::Ones();
	dcomplex k = 1.0;

	iso_bem_t iso_b(mesh);
	auto result1 = iso_b.eval(x0, k);
	std::cout << "isoparametric result : " << result1 << std::endl;

	const_bem_t const_b(mesh);
	auto result2 = const_b.eval(x0, k);
	std::cout << "constant result : " << result2 << std::endl;

	return 0;
}

