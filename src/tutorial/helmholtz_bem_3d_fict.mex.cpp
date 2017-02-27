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
#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
typedef NiHu::mex::complex_matrix<double> cMatrix;
//![Header]

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
//![Mesh]
	dMatrix
 		surf_nodes(rhs[0]), surf_elem(rhs[1]),	
		chief_nodes(rhs[2]), chief_elem(rhs[3]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::tria_1_tag(), NiHu::quad_1_tag());
	auto chief_mesh = NiHu::create_mesh(chief_nodes, chief_elem, NiHu::quad_1_tag());
//![Mesh]
	
//! [Function spaces]
	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &chief_sp = NiHu::dirac(NiHu::constant_view(chief_mesh));
//! [Function spaces]

//! [Matrices]
	int n = surf_sp.get_num_dofs();
	int m = chief_sp.get_num_dofs();
	cMatrix
		L_surf(n, n, lhs[0]), M_surf(n, n, lhs[1]),
		L_chief(m, n, lhs[2]), M_chief(m, n, lhs[3]),
		Mt_surf(n, n, lhs[4]), N_surf(n, n, lhs[5]);
//! [Matrices]
	
//! [Integral operators]
	double k = *mxGetPr(rhs[4]);
	auto I = NiHu::identity_integral_operator();
	auto L = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(k));
	auto M = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));
	auto Mt = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLPt_kernel<double>(k));
	auto N = NiHu::create_integral_operator(NiHu::helmholtz_3d_HSP_kernel<double>(k));
//! [Integral operators]

//! [System matrices]
	// conventional equations
	L_surf << dirac(surf_sp) * L[surf_sp];
	M_surf << dirac(surf_sp) * M[surf_sp]  +  dirac(surf_sp) * (-.5*I)[surf_sp];
	// CHIEF equations
	L_chief << chief_sp * L[surf_sp];
	M_chief << chief_sp * M[surf_sp];
	// Burton-Miller equations
	Mt_surf  << dirac(surf_sp) * Mt[surf_sp] +  dirac(surf_sp) * (.5*I)[surf_sp];
	N_surf  << dirac(surf_sp) * N[surf_sp];
//! [System matrices]
}

