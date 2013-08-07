#include "util/mex_matrix.hpp"
#include "bem/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 4 || nrhs < 4)
		throw("Too few input or output arguments");

	dMatrix surf_nodes(rhs[0]), surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag(), _tria_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	dMatrix field_nodes(rhs[2]), field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	auto L = create_integral_operator(helmholtz_3d_SLP_kernel(1.0));
	auto M = create_integral_operator(helmholtz_3d_DLP_kernel(1.0));
	auto I = -.5 * identity_integral_operator();

	auto n = surf_sp.get_num_dofs();
	cMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]);

	Ls << ( dirac(surf_sp) * L[surf_sp] );
	Ms << ( dirac(surf_sp) * M[surf_sp] ) + ( dirac(surf_sp) * I[surf_sp] );

	auto m = field_sp.get_num_dofs();
	cMatrix Lf(m, n, lhs[2]), Mf(m, n, lhs[3]);

	Lf << ( field_sp * L[surf_sp] );
	Mf << ( field_sp * M[surf_sp] );
}

