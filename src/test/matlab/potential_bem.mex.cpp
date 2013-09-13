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

#include "util/mex_matrix.hpp"
#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

#include <chrono>

typedef mex::real_matrix<double> dMatrix;


void mexFunction(int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 8 || nrhs < 4)
		throw("Too few input or output arguments");

	// create used integral operators

	auto L = create_integral_operator(laplace_3d_SLP_kernel());
	auto M = create_integral_operator(laplace_3d_DLP_kernel());
	auto Mt = create_integral_operator(laplace_3d_DLPt_kernel());
	auto N = create_integral_operator(laplace_3d_HSP_kernel());
	auto I = identity_integral_operator();

	// create surface mesh and test and trial function spaces

	dMatrix surf_nodes(rhs[0]), surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag(), _tria_1_tag());
	auto const &trial_sp = constant_view(surf_mesh);	// constant
	auto const &test_sp = dirac(trial_sp);				// collocation

	auto n = trial_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]);

	Ls << ( test_sp * L[trial_sp] );
	Ms << ( test_sp * M[trial_sp] ) + ( test_sp * (-.5*I)[trial_sp] );

	// create field point mesh

	dMatrix field_nodes(rhs[2]), field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[2]), Mf(m, n, lhs[3]), Mtf(m, n, lhs[4]), Nf(m, n, lhs[5]);

	// evaluate radiation matrices with separate kernel evaluation

	auto start = std::chrono::steady_clock::now();	// start timer
	Lf  << ( field_sp * L [trial_sp] );
	Mf  << ( field_sp * M [trial_sp] );
	Mtf << ( field_sp * Mt[trial_sp] );
	Nf  << ( field_sp * N [trial_sp] );
	auto stop = std::chrono::steady_clock::now();	// stop timer
	dMatrix dur_separate(1, 1, lhs[6]);
	dur_separate(0,0) = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

	// evaluate radiation matrices again with couple kernel evaluation

	auto CoupleOp = create_integral_operator(
		laplace_3d_SLP_kernel(),
		laplace_3d_DLP_kernel(),
		laplace_3d_DLPt_kernel(),
		laplace_3d_HSP_kernel()
	);
	start = std::chrono::steady_clock::now();	// start timer
	create_couple(Lf, Mf, Mtf, Nf)  << ( field_sp * CoupleOp [trial_sp] );
	stop = std::chrono::steady_clock::now();	// stop timer
	dMatrix dur_couple(1, 1, lhs[7]);
	dur_couple(0,0) = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
}

