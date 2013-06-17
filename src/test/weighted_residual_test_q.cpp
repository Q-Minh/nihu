#include "../bem/weighted_residual.hpp"

typedef tmp::vector<tria_1_elem> elem_vector;
typedef Mesh<elem_vector> mesh_t;
typedef function_space<mesh_t, constant_field, function_field> trial_space_t;
typedef trial_space_t test_space_t;	// Galerkin
typedef helmholtz_HG_kernel kernel_t;
typedef weighted_residual<kernel_t, test_space_t, trial_space_t> wr_t;

int main(void)
{
	try
	{
		Eigen::Matrix<double, Eigen::Dynamic, 3> nodes(3,3);
		nodes <<
			0.0, 0.0, 0.0,
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0;

		Eigen::Matrix<unsigned, Eigen::Dynamic, 3+1> elements(1,3+1);
		elements <<
			231, 0, 1, 2;

		mesh_t mesh(nodes, elements);

		test_space_t test_func(mesh);
		trial_space_t trial_func(mesh);
		wr_t wr(test_func, trial_func);

		kernel_t::set_wave_number(1.0);

		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> big_mat_t;
		big_mat_t a(test_func.get_num_dofs(), trial_func.get_num_dofs());
		a.setZero();

		big_mat_t b = a;

		couple<big_mat_t, big_mat_t> result(a, b);
		wr.eval(result);

		std::cout << result.first() << std::endl << std::endl;
		std::cout << result.second() << std::endl << std::endl;
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

