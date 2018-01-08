#include "util/mex_matrix.hpp"
#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::line_1_tag());
	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	int n = surf_sp.get_num_dofs();
	dMatrix L_surf(n, n, lhs[0]);
	dMatrix M_surf(n, n, lhs[1]);
	auto L = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto M = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());
	L_surf << dirac(surf_sp) * L[surf_sp];
	M_surf << dirac(surf_sp) * M[surf_sp];
}
