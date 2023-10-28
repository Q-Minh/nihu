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

typedef NiHu::field_dimension::_2d field_dim_t;
static const unsigned int field_dim = field_dim_t::value;
typedef NiHu::space_3d<> space_t;
static const unsigned int space_dim = space_t::dimension;

typedef NiHu::mex::real_matrix<double> dMatrix;

typedef Eigen::Matrix<double, field_dim, field_dim> field_var_t;
typedef Eigen::Matrix<double, space_dim, space_dim> space_var_t;

typedef NiHu::gaussian_covariance_kernel<space_t, field_dim_t> kernel_t;

// [D, B] = mex(nodes, elements, field_var, space_var);
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix nodes(rhs[0]), elem(rhs[1]);
	auto mesh = NiHu::create_mesh(nodes, elem, NiHu::tria_1_tag());
	auto const &w = NiHu::isoparametric_view(mesh, field_dim_t());

	
	field_var_t field_var = dMatrix(rhs[2]);
	space_var_t space_var = dMatrix(rhs[3]);
	kernel_t kernel(field_var, space_var);
	auto C = NiHu::create_integral_operator(kernel);
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

