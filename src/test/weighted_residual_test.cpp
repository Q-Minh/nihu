#include "../bem/weighted_residual.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;
typedef Mesh<elem_vector> mesh_t;
typedef function_space<mesh_t, constant_field, dirac_field> test_space_t;
typedef function_space<mesh_t, constant_field, function_field> trial_space_t;

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

	mesh_t mesh;
	for (unsigned i = 0; i < 6; ++i)
		mesh.add_node(c+i*3);
	for (unsigned i = 0; i < 3; ++i)
		mesh.add_elem(e+i*5);

	test_space_t test_func(mesh);
	trial_space_t trial_func(mesh);
	weighted_residual<helmholtz_G_kernel, test_space_t, trial_space_t> wr(test_func, trial_func);

	Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> a(test_func.get_num_dofs(), trial_func.get_num_dofs());
	a.setZero();

	wr.eval(a);

	std::cout << a << std::endl;

	return 0;
}

