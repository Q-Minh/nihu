#include "util/mex_matrix.hpp"
#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]), field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::line_1_tag());
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem, NiHu::line_1_tag());
	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &field_sp = NiHu::constant_view(field_mesh);
	size_t n = surf_sp.get_num_dofs();
	dMatrix L_surf(n, n, lhs[0]);
	dMatrix M_surf(n, n, lhs[1]);
	size_t m = field_sp.get_num_dofs();
	dMatrix L_field(m, n, lhs[2]);
	dMatrix M_field(m, n, lhs[3]);
	auto L = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto M = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());
	L_surf << dirac(surf_sp) * L[surf_sp];
	M_surf << dirac(surf_sp) * M[surf_sp];
	L_field << dirac(field_sp) * L[surf_sp];
	M_field << dirac(field_sp) * M[surf_sp];
}
