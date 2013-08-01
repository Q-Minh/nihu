#include "util/mex_matrix.hpp"
#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

typedef mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 4 || nrhs < 4)
		throw("Too few input or output arguments");

	dMatrix surf_nodes(rhs[0]), surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	dMatrix field_nodes(rhs[2]), field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	auto L = create_integral_operator(poisson_SLP_kernel());
	auto M = create_integral_operator(poisson_DLP_kernel());
	auto I = -.5 * identity_integral_operator();

	auto n = surf_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]);

	Ls << ( surf_sp * L[surf_sp] );
	Ms << ( surf_sp * M[surf_sp] ) + ( surf_sp * I[surf_sp] );

	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[2]), Mf(m, n, lhs[3]);

	Lf << ( field_sp * L[surf_sp] );
	Mf << ( field_sp * M[surf_sp] );
}

