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
#include "library/lib_element.hpp"
#include "library/elastostatics_kernel.hpp"
#include "library/elastostatics_singular_integrals.hpp"
#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
//![Header]

//![Mesh]
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::quad_2_tag());
//![Mesh]

//! [Function spaces]
	auto const &surf_sp = NiHu::constant_view(surf_mesh, NiHu::_3d());
//! [Function spaces]

//! [Matrices]
	int n = surf_sp.get_num_dofs();
	dMatrix Us(n, n, lhs[0]), Ts(n, n, lhs[1]);
//! [Matrices]

//! [Integral operators]
	double nu = *mxGetPr(rhs[2]);
	double mu = *mxGetPr(rhs[3]);
	auto U = NiHu::create_integral_operator(NiHu::elastostatics_3d_U_kernel(nu, mu));
	auto T = NiHu::create_integral_operator(NiHu::elastostatics_3d_T_kernel(nu, mu));
//! [Integral operators]

//! [System matrices]
	Us << NiHu::dirac(surf_sp) * U[surf_sp];
	Ts << NiHu::dirac(surf_sp) * T[surf_sp];
}
//! [System matrices]

