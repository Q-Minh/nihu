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

#include <boost/math/constants/constants.hpp>

#include "nihu/core/weighted_residual.hpp"
#include "nihu/library/elastodynamics_kernel.hpp"
#include "nihu/library/elastodynamics_singular_integrals.hpp"
#include "nihu/library/lib_element.hpp"
#include "interface/read_off_mesh.hpp"

#include <fstream>
#include <sstream>

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

int main(int argc, char *argv[])
{
	using namespace boost::math::double_constants;

	auto mesh = NiHu::read_off_mesh(argv[1], NiHu::tria_1_tag(), NiHu::quad_1_tag());
	auto const &v = NiHu::constant_view(mesh, NiHu::field_dimension::_3d());

	auto field = NiHu::read_off_mesh(argv[2], NiHu::quad_1_tag());
	auto const &w = NiHu::constant_view(field, NiHu::field_dimension::_3d());

	double nu = .33;
	double rho = 100;
	double mu = 1e8;
	double omega = two_pi*50.0;
	auto U_op = NiHu::create_integral_operator(NiHu::elastodynamics_3d_U_kernel(nu, rho, mu, omega));
	auto T_op = NiHu::create_integral_operator(NiHu::elastodynamics_3d_T_kernel(nu, rho, mu, omega));

	cMatrix Us(v.get_num_dofs(), v.get_num_dofs());
	Us.setZero();
	Us << dirac(v) * U_op[v];

	cMatrix Ts(v.get_num_dofs(), v.get_num_dofs());
	Ts.setZero();
	Ts << dirac(v) * T_op[v];

	cMatrix Uf(w.get_num_dofs(), v.get_num_dofs());
	Uf.setZero();
	Uf << dirac(w) * U_op[v];

	cMatrix Tf(w.get_num_dofs(), v.get_num_dofs());
	Tf.setZero();
	Tf << dirac(w) * T_op[v];

	std::cout << "ready" << std::endl;

	std::ofstream o;
	o.open("Us.mtx", std::ofstream::out); o << Us; o.close();
	o.open("Uf.mtx", std::ofstream::out); o << Uf; o.close();
	o.open("Ts.mtx", std::ofstream::out); o << Ts; o.close();
	o.open("Tf.mtx", std::ofstream::out); o << Tf; o.close();

	return 0;
}

