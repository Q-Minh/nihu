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

#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/lib_element.hpp"
#include <chrono>

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
//! [Straightforward]
	// import Matlab matrices
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]), field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elem, quad_1_tag());
	auto field_mesh = create_mesh(field_nodes, field_elem, quad_1_tag());

	// build function spaces
	auto const &surf_sp = isoparametric_view(surf_mesh);
	auto const &field_sp = dirac(constant_view(field_mesh));

	// number of DOF and radiation points
	int n = surf_sp.get_num_dofs();
	int m = field_sp.get_num_dofs();
//! [Straightforward]

//! [Without coupling]
	// preallocation
	cMatrix Lmat(m, n, lhs[0]), Mmat(m, n, lhs[1]), Mtmat(m, n, lhs[2]), Nmat(m, n, lhs[3]);

	double k = *mxGetPr(rhs[4]);	// get wave number from Matab

	// integral operators
	auto L = create_integral_operator(helmholtz_3d_SLP_kernel<double>(k));
	auto M = create_integral_operator(helmholtz_3d_DLP_kernel<double>(k));
	auto Mt = create_integral_operator(helmholtz_3d_DLPt_kernel<double>(k));
	auto N = create_integral_operator(helmholtz_3d_HSP_kernel<double>(k));

	auto start = std::chrono::steady_clock::now();	// start timer

	Lmat  << field_sp * L[surf_sp];
	Mmat  << field_sp * M[surf_sp];
	Mtmat << field_sp * Mt[surf_sp];
	Nmat  << field_sp * N[surf_sp];

	auto stop = std::chrono::steady_clock::now();	// stop timer
	dMatrix dur_plain(1, 1, lhs[8]);				// pass duration to Matlab
	dur_plain(0, 0) = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
//! [Without coupling]

//! [With coupling]
	// preallocation

	cMatrix Lmat_c(m, n, lhs[4]), Mmat_c(m, n, lhs[5]), Mtmat_c(m, n, lhs[6]), Nmat_c(m, n, lhs[7]);

	start = std::chrono::steady_clock::now();	// start timer

	// create couple operator
	auto op = create_integral_operator(
		helmholtz_3d_SLP_kernel<double>(k),
		helmholtz_3d_DLP_kernel<double>(k),
		helmholtz_3d_DLPt_kernel<double>(k),
		helmholtz_3d_HSP_kernel<double>(k));

	// evaluate into couple_matrix
	create_couple(Lmat_c, Mmat_c, Mtmat_c, Nmat_c) << field_sp * op[surf_sp];

	stop = std::chrono::steady_clock::now();	// stop timer
	dMatrix dur_coup(1, 1, lhs[9]);				// pass duration to Matlab
	dur_coup(0, 0) = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
//! [With coupling]
}

