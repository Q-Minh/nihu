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

#include "../core/integral_operator.hpp"
#include "../tmp/vector.hpp"
#include "../core/mesh.hpp"
#include "../core/function_space.hpp"
#include "../core/projection.hpp"
#include "../core/weighted_residual.hpp"

#include "../library/unit_kernel.hpp"
#include "../library/poisson_kernel.hpp"
#include "../library/helmholtz_kernel.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(9,3);
	nodes <<
		-1.0, -1.0, 0.0,
		 0.0, -1.0, 0.0,
		 1.0, -1.0, 0.0,

		-1.0,  0.0, 0.0,
		 0.0,  0.0, 0.0,
		 1.0,  0.0, 0.0,

		-1.0,  1.0, 0.0,
		 0.0,  1.0, 0.0,
		 1.0,  1.0, 0.0;

	uMatrix elements(5, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	mesh_t msh(nodes, elements);

	auto trial = constant_view(msh);	// constant
	auto const &test = dirac(trial);
	int nDOF = trial.get_num_dofs();
	dMatrix A(nDOF, nDOF);
	A.setZero();

	auto b_op = create_integral_operator(poisson_G_kernel());
	auto id_op = .5 * identity_integral_operator();
	auto proj = b_op[trial] + id_op[trial];
	auto wr = test * proj;
	wr.eval(A);

	std::cout << A << std::endl;

	return 0;
}

