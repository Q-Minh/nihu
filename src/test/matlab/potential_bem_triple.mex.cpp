#include "util/mex_matrix.hpp"
#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

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
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	dMatrix field_nodes(rhs[2]);
	dMatrix field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	// generating integral operators

	auto LM = create_integral_operator(
		create_couple_kernel(
			poisson_3d_SLP_kernel(),
			poisson_3d_DLP_kernel(),
			poisson_3d_DLP_kernel())
		);
	auto I = -.5 * identity_integral_operator();

	// surface system matrices
	
	auto n = surf_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]);
	dMatrix Ms(n, n, lhs[1]);
	dMatrix Ms2(n, n, lhs[2]);

	( surf_sp * LM[surf_sp] ).eval( create_couple(Ls, Ms, Ms2) );
	( surf_sp *  I[surf_sp] ).eval( Ms );
	( surf_sp *  I[surf_sp] ).eval( Ms2 );

	// field point system matrices
	
	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[3]);
	dMatrix Mf(m, n, lhs[4]);
	dMatrix Mf2(m, n, lhs[5]);

	( field_sp * LM[surf_sp] ).eval( create_couple(Lf, Mf, Mf2) );
}

