/* $Make: mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" mass_matrix.mex.cpp -I../../ -I../../../../eigen -output mass_matrix $ */

#include "util/mex_matrix.h"
#include "bem/weighted_residual.hpp"

void mexFunction(
	int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	mex::real_matrix nodes(rhs[0]);
	mex::real_matrix elements(rhs[1]);

	mesh<tmp::vector<quad_1_elem> > msh(nodes, elements);
	auto f_sp = isoparametric_view(msh);
	unsigned ndof = f_sp.get_num_dofs();
	mex::real_matrix result(ndof, ndof, lhs[0]);
	auto I = identity_integral_operator();

	( f_sp * I[f_sp] ).eval(result);
}

