#include "../bem/weighted_residual.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;
typedef Mesh<elem_vector> mesh_t;
typedef function_space<mesh_t, constant_field, dirac_field> test_space_t;
typedef function_space<mesh_t, constant_field, function_field> trial_space_t;
typedef helmholtz_HG_kernel kernel_t;
typedef weighted_residual<kernel_t, test_space_t, trial_space_t> wr_t;

int main(void)
{
	try
	{
		int const n = 50;

		Eigen::Matrix<double, Eigen::Dynamic, 3> nodes((n+1)*(n+1),3);
		Eigen::Matrix<unsigned, Eigen::Dynamic, 5> elements(n*n,5);

		for (int i = 0; i < (n+1); ++i)
		{
			for (int j = 0; j < (n+1); ++j)
			{
				int row = i*(n+1)+j;
				nodes(row,0) = double(i);
				nodes(row,1) = double(j);
				nodes(row,2) = 0.0;
			}
		}

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				int row = i*n+j;
				elements(row,0) = 4;
				elements(row,1) = i*(n+1)+j+0;
				elements(row,2) = i*(n+1)+j+1;
				elements(row,3) = i*(n+1)+j+n+2;
				elements(row,4) = i*(n+1)+j+n+1;
			}
		}

		mesh_t mesh;
		mesh.build_from_mex<5>(nodes, elements);

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

	/*
		std::cout << result.first() << std::endl;
		std::cout << result.second() << std::endl;
		*/
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

