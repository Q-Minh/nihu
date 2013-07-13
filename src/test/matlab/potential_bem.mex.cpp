/* $Make: mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" potential_bem.mex.cpp -I../../ -I../../../../eigen -output potential_bem $ */

#include "util/mex_matrix.h"
#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

void mexFunction(
	int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 4 || nrhs < 4)
		return;
	mex::real_matrix nodes(rhs[0]);
	mex::real_matrix elements(rhs[1]);

	mex::real_matrix field_nodes(rhs[2]);
	mex::real_matrix field_elements(rhs[3]);

	mesh<tmp::vector<quad_1_elem> > surf_mesh(nodes, elements);
	mesh<tmp::vector<quad_1_elem> > field_mesh(field_nodes, field_elements);

	auto const &surf_sp = isoparametric_view(surf_mesh);
	auto const &field_sp = dirac(constant_view(field_mesh));

	unsigned N = surf_sp.get_num_dofs();
	unsigned M = field_sp.get_num_dofs();

	mex::real_matrix Gs(N, N, lhs[0]);
	mex::real_matrix Hs(N, N, lhs[1]);

	mex::real_matrix Gf(M, N, lhs[2]);
	mex::real_matrix Hf(M, N, lhs[3]);

	auto I = identity_integral_operator();
	auto G = create_integral_operator(poisson_G_kernel());
	auto H = create_integral_operator(poisson_H_kernel());

	( surf_sp * G[surf_sp] ).eval(Gs);
	( surf_sp * (H[surf_sp] + (-.5*I)[surf_sp]) ).eval(Hs);

	( field_sp * G[surf_sp] ).eval(Gf);
	( field_sp * H[surf_sp] ).eval(Hf);
}

