#include "util/mex_matrix.hpp"
#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elem, tria_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	dMatrix field_nodes(rhs[2]), field_elem(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elem, tria_1_tag());
	auto const &field_sp = constant_view(field_mesh);

	int n = surf_sp.get_num_dofs();
	dMatrix
		L_surf(n, n, lhs[0]), M_surf(n, n, lhs[1]);
	int m = field_sp.get_num_dofs();
	dMatrix
		L_field(m, n, lhs[2]), M_field(m, n, lhs[3]);

	auto L = create_integral_operator(laplace_3d_SLP_kernel());
	auto M = create_integral_operator(laplace_3d_DLP_kernel());

	L_surf << dirac(surf_sp) * L[surf_sp];
	M_surf << dirac(surf_sp) * M[surf_sp];
	L_field << dirac(field_sp) * L[surf_sp];
	M_field << dirac(field_sp) * M[surf_sp];
}

