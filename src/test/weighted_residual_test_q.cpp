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

typedef field<quad_1_elem, quad_0_shape_set> quad_1_field;

typedef tmp::vector<quad_1_field> field_vector;
typedef function_space<field_vector> function_space_t;
typedef unit_kernel kernel_t;
typedef weighted_residual<true, kernel_t, function_space_t, function_space_t> wr_t;

int main(void)
{
	try
	{
		Eigen::Matrix<double, Eigen::Dynamic, 3> nodes(10,3);
		nodes <<
			0.0, 0.0, 0.0,
			1.0, 0.0, 0.0,
			2.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			1.0, 1.0, 0.0,
			2.0, 1.0, 0.0,
			0.0, 2.0, 0.0,
			1.0, 2.0, 0.0,
			2.0, 2.0, 0.0,
			3.0, 2.0, 0.0;

		Eigen::Matrix<unsigned, Eigen::Dynamic, 4+1+1> fields(3,4+1+1);
		fields <<
			quad_1_field::id, 0, 1, 4, 3, 0,
			quad_1_field::id, 3, 4, 7, 6, 1,
			quad_1_field::id, 4, 5, 8, 7, 2;

		function_space_t fsp(nodes, fields);

		// kernel_t::set_wave_number(1.0);

		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> big_mat_t;
		big_mat_t a(fsp.get_num_dofs(), fsp.get_num_dofs());
		a.setZero();

		//		big_mat_t b = a;

		//		couple<big_mat_t, big_mat_t> result(a, b);
		wr_t::eval(a, fsp, fsp);

		std::cout << a << std::endl << std::endl;
		//		std::cout << result.second() << std::endl << std::endl;
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

