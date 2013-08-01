/* $Make: mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" potential_bem.mex.cpp -I../../ -I/usr/local/include/eigen3 -output potential_bem $ */

#include "util/mex_matrix.h"
#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

typedef mex::real_matrix<double> dMatrix;

void mexFunction(
	int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 4 || nrhs < 4)
		return;

	// generating function spaces

	dMatrix surf_nodes(rhs[0]);
	dMatrix surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	dMatrix field_nodes(rhs[2]);
	dMatrix field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	// generating integral operators

	auto LM = create_integral_operator(
		create_couple_kernel(
			poisson_SLP_kernel(),
			poisson_DLP_kernel()));
	auto I = -.5 * identity_integral_operator();

	// surface system matrices

	auto n = surf_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]);
	dMatrix Ms(n, n, lhs[1]);

	(Ls, Ms) << ( surf_sp * LM[surf_sp] );
	Ms << ( surf_sp *  I[surf_sp] );

	// field point system matrices

	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[2]);
	dMatrix Mf(m, n, lhs[3]);

	(Lf, Mf) << ( field_sp * LM[surf_sp] );
}

