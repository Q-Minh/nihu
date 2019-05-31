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
#include "library/covariance_kernel.hpp"
#include "util/mex_matrix.hpp"
#include "library/lib_element.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

// [D, B] = mex(nodes, elements, sigma, d);
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix nodes(rhs[0]), elem(rhs[1]);
	auto mesh = NiHu::create_mesh(nodes, elem, NiHu::tria_1_tag());
	auto const &w = NiHu::isoparametric_view(mesh);

	double sigma = NiHu::mex::get_scalar<double>(rhs[2]);
	double d = NiHu::mex::get_scalar<double>(rhs[3]);
	auto C = NiHu::create_integral_operator(NiHu::covariance_kernel<NiHu::space_3d<> >(sigma, d));
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

