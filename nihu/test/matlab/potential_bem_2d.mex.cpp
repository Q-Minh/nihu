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

#include "nihu/util/mex_matrix.hpp"
#include "nihu/core/weighted_residual.hpp"
#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/laplace_nearly_singular_integrals.hpp"
#include "nihu/library/laplace_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 6 || nrhs < 4)
		throw("Too few input or output arguments");

	dMatrix surf_nodes(rhs[0]), surf_elements(rhs[1]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elements, NiHu::line_1_tag());
	auto const &trial_sp = NiHu::constant_view(surf_mesh);
	auto const &test_sp = NiHu::dirac(trial_sp);

	dMatrix field_nodes(rhs[2]), field_elements(rhs[3]);
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elements, NiHu::line_1_tag());
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field_mesh));

	auto L = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto M = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());
	auto Mt = NiHu::create_integral_operator(NiHu::laplace_2d_DLPt_kernel());
	auto N = NiHu::create_integral_operator(NiHu::laplace_2d_HSP_kernel());
	auto I = NiHu::identity_integral_operator();

	auto n = trial_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]), Mts(n, n, lhs[2]), Ns(n, n, lhs[3]);

	Ls << ( test_sp * L[trial_sp] );
	Ms << ( test_sp * M[trial_sp] ) + ( test_sp * (-.5*I)[trial_sp] );
	Mts << ( test_sp * Mt[trial_sp] ) + ( test_sp * (.5*I)[trial_sp] );
	Ns << ( test_sp * N[trial_sp] );

	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[4]), Mf(m, n, lhs[5]);

	Lf << (field_sp * L[trial_sp]);
	Mf << (field_sp * M[trial_sp]);
}

