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

#include "../bem/integral_operator.hpp"
#include "../library/unit_kernel.hpp"

typedef tmp::vector<line_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(3,2);
	nodes <<
		-1.0, 0.0,
		 0.0, 0.0,
		 1.0, 0.0;

	uMatrix elements(2, 1+2);
	elements <<
		line_1_elem::id, 0, 1,
		line_1_elem::id, 1, 2;

	mesh_t msh(nodes, elements);

	auto f_sp = create_function_space_view(msh, field_option::constant());

	dMatrix A(f_sp.get_num_dofs(), f_sp.get_num_dofs());
	A.setZero();
	
	auto b_op = create_integral_operator(unit_kernel<space_2d>(), operator_option::local());
	A += b_op(f_sp, f_sp, formalism::full_dirac());

	std::cout << A << std::endl << std::endl;
	std::cout << b_op.get_kernel().get_num_evaluations() << std::endl;

	return 0;
}

