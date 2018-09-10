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
#include "util/mex_matrix.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
typedef NiHu::mex::complex_matrix<double> cMatrix;

// [Ls, Ms, Mts, Ns, Lf, Mf, Mtf, Nf] = mex(s_nodes, s_elements, f_nodes, f_elements, k);
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	dMatrix field_nodes(rhs[2]), field_elem(rhs[3]);

	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::tria_1_tag(), NiHu::quad_1_tag());
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem, NiHu::tria_1_tag(), NiHu::quad_1_tag());
	
	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &field_sp = NiHu::constant_view(field_mesh);

	int n = surf_sp.get_num_dofs();
	int m = field_sp.get_num_dofs();
	
	cMatrix
		L_surf(n, n, lhs[0]), M_surf(n, n, lhs[1]),
		Mt_surf(n, n, lhs[2]), N_surf(n, n, lhs[3]),
		L_field(m, n, lhs[4]), M_field(m, n, lhs[5]),
		Mt_field(m, n, lhs[6]), N_field(m, n, lhs[7]);

	// Retrieve wave number
	double k = *mxGetPr(rhs[4]);
	auto I = NiHu::identity_integral_operator();
	auto L = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(k));
	auto M = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));
	auto Mt = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLPt_kernel<double>(k));
	auto N = NiHu::create_integral_operator(NiHu::helmholtz_3d_HSP_kernel<double>(k));

	// conventional equations
	L_surf << NiHu::dirac(surf_sp) * L[surf_sp];
	M_surf << NiHu::dirac(surf_sp) * M[surf_sp]
		+  NiHu::dirac(surf_sp) * (-.5*I)[surf_sp];

	// Burton-Miller equations
	Mt_surf  << NiHu::dirac(surf_sp) * Mt[surf_sp]
		+ NiHu::dirac(surf_sp) * (.5*I)[surf_sp];
	N_surf  << NiHu::dirac(surf_sp) * N[surf_sp];
	
	// Field matrices
	L_field << NiHu::dirac(field_sp) * L[surf_sp];
	M_field << NiHu::dirac(field_sp) * M[surf_sp];
	
	// normal derived Field matrices
	Mt_field << NiHu::dirac(field_sp) * Mt[surf_sp];
	N_field << NiHu::dirac(field_sp) * N[surf_sp];
}
