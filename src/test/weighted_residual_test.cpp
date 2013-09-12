// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
		int const n = 2;

		Eigen::Matrix<double, Eigen::Dynamic, 3> nodes((n+1)*(n+1),3);
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

		Eigen::Matrix<unsigned, Eigen::Dynamic, 5> elements(n*n,5);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				int row = i*n+j;
				elements(row,0) = 241;
				elements(row,1) = i*(n+1)+j+0;
				elements(row,2) = i*(n+1)+j+1;
				elements(row,3) = i*(n+1)+j+n+2;
				elements(row,4) = i*(n+1)+j+n+1;
			}
		}

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

