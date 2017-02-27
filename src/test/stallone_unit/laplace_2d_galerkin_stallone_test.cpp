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
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(4,2);
	uMatrix elements(4,3);

	nodes <<
		0.0, 0.0,
		1.0, 0.0,
		1.0, 1.0,
		0.0, 1.0;
	elements <<
		NiHu::line_1_elem::id, 0, 1,
		NiHu::line_1_elem::id, 1, 2,
		NiHu::line_1_elem::id, 2, 3,
		NiHu::line_1_elem::id, 3, 0;

	auto mesh = NiHu::create_mesh(nodes, elements, NiHu::line_1_tag());
	auto const &fspace = NiHu::constant_view(mesh, NiHu::_1d());
	auto nDof = fspace.get_num_dofs();

	dMatrix L(nDof, nDof), M(nDof, nDof);
	L.setZero(); M.setZero();

	auto L_op = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto M_op = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());

	L << fspace * L_op[fspace];
	M << fspace * M_op[fspace];

	std::cout << L << std::endl;
	std::cout << M << std::endl;

	return 0;
}
