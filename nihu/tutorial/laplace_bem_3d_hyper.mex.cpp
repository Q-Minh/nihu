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

//![Header]
#include "nihu/util/mex_matrix.hpp"
#include "nihu/core/weighted_residual.hpp"
#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/laplace_nearly_singular_integrals.hpp"
#include "nihu/library/laplace_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
//![Header]

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
//![Mesh]
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::tria_1_tag(), NiHu::quad_1_tag());
//![Mesh]

//! [Function spaces]
	auto const &surf_sp = NiHu::constant_view(surf_mesh);
//! [Function spaces]

//! [Matrices]
	size_t n = surf_sp.get_num_dofs();
	dMatrix
		L_surf(n, n, lhs[0]), M_surf(n, n, lhs[1]),
		Mt_surf(n, n, lhs[2]), N_surf(n, n, lhs[3]);
//! [Matrices]

//! [Integral operators]
	auto I = NiHu::identity_integral_operator();
	auto L = NiHu::create_integral_operator(NiHu::laplace_3d_SLP_kernel());
	auto M = NiHu::create_integral_operator(NiHu::laplace_3d_DLP_kernel());
	auto Mt = NiHu::create_integral_operator(NiHu::laplace_3d_DLPt_kernel());
	auto N = NiHu::create_integral_operator(NiHu::laplace_3d_HSP_kernel());
//! [Integral operators]

//! [System matrices]
	// conventional equations
	L_surf << dirac(surf_sp) * L[surf_sp];
	M_surf << dirac(surf_sp) * M[surf_sp]  +  dirac(surf_sp) * (-.5*I)[surf_sp];
	// hypersingular equations
	Mt_surf  << dirac(surf_sp) * Mt[surf_sp] +  dirac(surf_sp) * (.5*I)[surf_sp];
	N_surf  << dirac(surf_sp) * N[surf_sp];
//! [System matrices]
}

