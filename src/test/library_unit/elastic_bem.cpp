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

#include <iostream>
#include "core/weighted_residual.hpp"
#include "library/lib_element.hpp"
#include "library/elastostatics_kernel.hpp"
#include "library/elastostatics_singular_integrals.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(4,3);
	nodes <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0,
		0.0, 1.0, 0.0;
	uMatrix elements(1, 5);
	elements << quad_1_elem::id, 0, 1, 2, 3;
	auto mesh = create_mesh(nodes, elements, quad_1_tag());
	
	nodes <<
		0.0, 0.0, 2.0,
		1.0, 0.0, 2.0,
		1.0, 1.0, 2.0,
		0.0, 1.0, 2.0;
	auto field = create_mesh(nodes, elements, quad_1_tag());
	
	auto K = create_integral_operator(elastostatics_3d_U_kernel(.33));
	
	auto const &w = isoparametric_view(field, _3d());
	auto const &v = constant_view(mesh, _3d());
	
	dMatrix res(w.get_num_dofs(), v.get_num_dofs());

	res.setZero();
	res << dirac(w) * K[v];
	std::cout << res << std::endl;
	
	return 0;
}
