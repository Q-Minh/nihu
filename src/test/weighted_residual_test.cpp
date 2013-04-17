#include "../bem/weighted_residual.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;
typedef Mesh<elem_vector> mesh_t;
typedef function_space<mesh_t, isoparametric_field, function_field> test_space_t;
typedef function_space<mesh_t, constant_field, function_field> trial_space_t;
typedef helmholtz_HG_kernel kernel_t;
typedef weighted_residual<kernel_t, test_space_t, trial_space_t> wr_t;

#include <iostream>

int main(void)
{
	try
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
		wr_t wr(test_func, trial_func);

		kernel_t::set_wave_number(1.0);

		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> a(test_func.get_num_dofs(), trial_func.get_num_dofs());
		a.setZero();

		std::cout << wr.eval(a) << std::endl;
		std::cout << kernel_t::get_num_evaluations() << std::endl;
	}
	catch(const char *e)
	{
		std::cerr << e << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unhandled exception occurred" << std::endl;
	}

	return 0;
}

