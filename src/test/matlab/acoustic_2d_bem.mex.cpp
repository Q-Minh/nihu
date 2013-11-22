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
#include "library/helmholtz_singular_integrals.hpp"

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elem, line_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);
	int n = surf_sp.get_num_dofs();

	dMatrix field_nodes(rhs[2]), field_elem(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elem, line_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));
	int m = field_sp.get_num_dofs();

	cMatrix Ls(n, n, lhs[0]), Lf(m, n, lhs[1]);

	double k = *mxGetPr(rhs[4]);
	auto L = create_integral_operator(helmholtz_2d_SLP_kernel<double>(k));

	Ls << dirac(surf_sp) * L[surf_sp];
	Lf  << field_sp * L[surf_sp];
}

