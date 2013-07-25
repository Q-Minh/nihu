/* $Make: mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" potential_bem.mex.cpp -I../../ -I/usr/local/include/eigen3 -output potential_bem $ */

#include "util/mex_matrix.h"
#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

void mexFunction(
	int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 4 || nrhs < 4)
		return;

	// generating function spaces

	mex::real_matrix surf_nodes(rhs[0]);
	mex::real_matrix surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	mex::real_matrix field_nodes(rhs[2]);
	mex::real_matrix field_elements(rhs[3]);
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
	mex::real_matrix Ls(n, n, lhs[0]);
	mex::real_matrix Ms(n, n, lhs[1]);

	( surf_sp * LM[surf_sp] ).eval( create_couple(Ls, Ms) );
	( surf_sp *  I[surf_sp] ).eval( Ms );

	// field point system matrices
	
	auto m = field_sp.get_num_dofs();
	mex::real_matrix Lf(m, n, lhs[2]);
	mex::real_matrix Mf(m, n, lhs[3]);

	( field_sp * LM[surf_sp] ).eval( create_couple(Lf, Mf) );
}

