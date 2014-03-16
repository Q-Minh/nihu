// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
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

#include "core/weighted_residual.hpp"
#include "library/elastostatics_kernel.hpp"
#include "library/elastostatics_singular_integrals.hpp"
#include "library/lib_element.hpp"
#include "interface/read_off_mesh.hpp"

#include <fstream>
#include <sstream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

int main(int argc, char *argv[])
{
	auto mesh = read_off_mesh(argv[1], tria_1_tag(), quad_1_tag());
	auto const &v = constant_view(mesh, _3d());

	auto field = read_off_mesh(argv[2], quad_1_tag());
	auto const &w = constant_view(field, _3d());

	double nu = .33;
	auto U_op = create_integral_operator(elastostatics_3d_U_kernel(nu));
	auto T_op = create_integral_operator(elastostatics_3d_T_kernel(nu));

	dMatrix Us(v.get_num_dofs(), v.get_num_dofs());
	Us.setZero();
	Us << dirac(v) * U_op[v];

	dMatrix Ts(v.get_num_dofs(), v.get_num_dofs());
	Ts.setZero();
	Ts << dirac(v) * T_op[v];

	dMatrix Uf(w.get_num_dofs(), v.get_num_dofs());
	Uf.setZero();
	Uf << dirac(w) * U_op[v];

	dMatrix Tf(w.get_num_dofs(), v.get_num_dofs());
	Tf.setZero();
	Tf << dirac(w) * T_op[v];

	std::cout << "ready" << std::endl;

	std::ofstream o1;

	o1.open("Us.mtx", std::ofstream::out); o1 << Us; o1.close();
	o1.open("Uf.mtx", std::ofstream::out); o1 << Uf; o1.close();
	o1.open("Ts.mtx", std::ofstream::out); o1 << Ts; o1.close();
	o1.open("Tf.mtx", std::ofstream::out); o1 << Tf; o1.close();

	return 0;
}

