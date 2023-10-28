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
#include "nihu/library/covariance_kernel.hpp"
#include "nihu/util/mex_matrix.hpp"
#include "nihu/library/lib_element.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

typedef NiHu::exponential_covariance_kernel<NiHu::space_2d<>, NiHu::field_dimension::_1d> kernel_t;

/** @todo write usage for mex */

// [D, B] = mex(nodes, elements, var, d);
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	mexPrintf("%d \n", NiHu::line_1_elem::id);
	dMatrix nodes(rhs[0]), elem(rhs[1]);
	auto mesh = NiHu::create_mesh(nodes, elem, NiHu::line_1_tag());
	auto const &w = NiHu::constant_view(mesh);

	kernel_t::field_variance_t var = dMatrix(rhs[2]);
	double d = NiHu::mex::get_scalar<double>(rhs[3]);
	auto C = NiHu::create_integral_operator(kernel_t(var, d));
	auto I = NiHu::identity_integral_operator();
	
	size_t N = w.get_num_dofs();
	mexPrintf("Number of DOFs: %d\n", N);

	dMatrix D(N, N, lhs[0]), B(N, N, lhs[1]);
	D.setZero();
	B.setZero();
	mexPrintf("Matrices initialised\n");
	D << (w * C[w]);
	mexPrintf("D matrix ready\n");
	B << (w * I[w]);
	mexPrintf("B matrix ready\n");
}

