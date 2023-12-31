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

#include "nihu/core/weighted_residual.hpp"
#include "nihu/util/mex_matrix.hpp"
#include "nihu/library/helmholtz_kernel.hpp"
#include "nihu/library/helmholtz_nearly_singular_integrals.hpp"
#include "nihu/library/helmholtz_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
typedef NiHu::mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::tria_1_tag());
	dMatrix field_nodes(rhs[2]), field_elem(rhs[3]);
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem, NiHu::quad_1_tag());
	
//! [Function spaces]
	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field_mesh));
//! [Function spaces]

	size_t n = surf_sp.get_num_dofs();
	size_t m = field_sp.get_num_dofs();
	cMatrix Ms(n, n, lhs[0]), Mf(m, n, lhs[1]);
	
//! [System]
	double k = *mxGetPr(rhs[4]);
	auto I = NiHu::identity_integral_operator();
	auto M = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));

	Ms << dirac(surf_sp) * M[surf_sp] + dirac(surf_sp) * (.5*I)[surf_sp];
	Mf << field_sp * M[surf_sp];
//! [System]
}

