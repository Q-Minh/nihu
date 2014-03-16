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

#include "library/matsumoto_2010.hpp"
//#include "library/helmholtz_singular_integrals.hpp"
#include "core/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"


// dynamically resizeable double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// dynamically resizeable complex matrix
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
// dynamically resizeable unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main()
{
	// nodal coordinates in 3D, 3 nodes
	dMatrix nodes(3, 3);
	nodes <<
		-1.0, -1.0, 0.0,
		0.0, -1.0, 0.0,
		0.0, 0.0, 0.0;

	// element nodal indices
	uMatrix elements(1, 1 + 4);
	elements <<
		tria_1_elem::id, 0, 1, 2, 0;

	// create the mesh
	auto msh = create_mesh(nodes, elements, tria_1_tag());

	// create a piecewise constant function space
	auto const &w = constant_view(msh);

	// compute number of DOF and allocate result matrix
	int nDOF = w.get_num_dofs();
	cMatrix H_mat(nDOF, nDOF);
	H_mat.setZero();

	cMatrix G_mat(nDOF, nDOF);
	G_mat.setZero();

	double k = 0.05;

	// create integral operator from kernel and perform weighted double integral
	auto H = create_integral_operator(helmholtz_3d_HSP_kernel<double>(k));
	H_mat << dirac(w) * H[w];

	auto G = create_integral_operator(helmholtz_3d_SLP_kernel<double>(k));
	G_mat << dirac(w) * G[w];

	std::cout << "H = " << std::endl << H_mat << std::endl << std::endl;
	std::cout << "G = " << std::endl << G_mat << std::endl << std::endl;

	return 0;
}
