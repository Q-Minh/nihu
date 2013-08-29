#include "util/mex_matrix.hpp"
#include "bem/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"

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
			laplace_3d_SLP_kernel(),
			laplace_3d_DLP_kernel());
	auto I = identity_integral_operator();

	// surface system matrices

	auto n = surf_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]);
	dMatrix Ms(n, n, lhs[1]);

	(Ls, Ms) << ( surf_sp * LM[surf_sp] );
	Ms << ( surf_sp * (-.5*I)[surf_sp] );

	// field point system matrices

	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[2]);
	dMatrix Mf(m, n, lhs[3]);

	(Lf, Mf) << ( field_sp * LM[surf_sp] );
}

