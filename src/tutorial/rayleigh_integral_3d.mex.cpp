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

//! [Header]
#include "core/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "util/mex_matrix.hpp"
#include "library/lib_element.hpp"

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;
//! [Header]

//! [Mex function]
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
//! [Mex function]
{
//! [Meshes]
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]), field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elem, quad_1_tag());
	auto field_mesh = create_mesh(field_nodes, field_elem, quad_1_tag());
//! [Meshes]

//! [Function spaces]
	auto const &w = isoparametric_view(surf_mesh);
	auto const &dirac_f = dirac(constant_view(field_mesh));
	auto const &dirac_s = dirac(constant_view(surf_mesh));
//! [Function spaces]

//! [Kernel]
	double wave_number = *mxGetPr(rhs[4]);
	auto G = create_integral_operator(helmholtz_3d_SLP_kernel<double>(wave_number));
//! [Kernel]

//! [Matrices]
	cMatrix Z_trans(dirac_f.get_num_dofs(), w.get_num_dofs(), lhs[0]);
	cMatrix Z_rad(dirac_s.get_num_dofs(), w.get_num_dofs(), lhs[1]);
//! [Matrices]

//! [Weighted residual]
	Z_trans << (dirac_f * G[w]);
	Z_rad << (dirac_s * G[w]);
//! [Weighted residual]
}

