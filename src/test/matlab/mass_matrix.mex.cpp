/* $Make: mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" mass_matrix.mex.cpp -I../../ -I../../../../eigen -output mass_matrix $ */

#include "util/mex_matrix.hpp"
#include "core/weighted_residual.hpp"

void mexFunction(
	int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 2 || nrhs < 1)
		return;
	mex::real_matrix<double> nodes(rhs[0]);
	mex::real_matrix<double> elements(rhs[1]);

	mesh<tmp::vector<quad_1_elem> > msh(nodes, elements);
	auto f_sp = isoparametric_view(msh);
	unsigned ndof = f_sp.get_num_dofs();
	mex::real_matrix<double> result(ndof, ndof, lhs[0]);
	auto I = identity_integral_operator();

	result << ( f_sp * I[f_sp] );
}

