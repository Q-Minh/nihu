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
#include "library/lib_element.hpp"


typedef mex::real_matrix<double> dMatrix;

void mexFunction(
	int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 6 || nrhs < 4)
		return;

	// generating function spaces

	dMatrix surf_nodes(rhs[0]);
	dMatrix surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, quad_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	dMatrix field_nodes(rhs[2]);
	dMatrix field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	// generating integral operators

	auto LM = create_integral_operator(
			laplace_3d_SLP_kernel(),
			laplace_3d_DLP_kernel(),
			laplace_3d_DLP_kernel()
		);
	auto I = identity_integral_operator();

	// surface system matrices
	
	auto n = surf_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]);
	dMatrix Ms(n, n, lhs[1]);
	dMatrix Ms2(n, n, lhs[2]);

	create_couple(Ls, Ms, Ms2) << ( surf_sp * LM[surf_sp] );
	Ms << ( surf_sp *  (-.5*I)[surf_sp] );
	Ms2 << ( surf_sp *  (-.5*I)[surf_sp] );

	// field point system matrices
	
	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[3]);
	dMatrix Mf(m, n, lhs[4]);
	dMatrix Mf2(m, n, lhs[5]);

	create_couple(Lf, Mf, Mf2) << ( field_sp * LM[surf_sp] );
}

