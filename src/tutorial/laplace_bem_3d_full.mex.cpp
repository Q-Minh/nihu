#include "util/mex_matrix.hpp"
#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]), field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::tria_1_tag());
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem, NiHu::quad_1_tag());

	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &field_sp = NiHu::constant_view(field_mesh);

	int n = surf_sp.get_num_dofs();
	int m = field_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]), Mts(n, n, lhs[2]), Ds(n, n, lhs[3]),
		Gxxs(n, n, lhs[4]), 
		Lf(m, n, lhs[5]), Mf(m, n, lhs[6]), Mtf(m, n, lhs[7]), Df(m, n, lhs[8]);

	auto L = NiHu::create_integral_operator(NiHu::laplace_3d_SLP_kernel());
	auto M = NiHu::create_integral_operator(NiHu::laplace_3d_DLP_kernel());
	auto Mt = NiHu::create_integral_operator(NiHu::laplace_3d_DLPt_kernel());
	auto D = NiHu::create_integral_operator(NiHu::laplace_3d_HSP_kernel());
	auto Gxx = NiHu::create_integral_operator(NiHu::laplace_3d_Gxx_kernel());

	Ls << NiHu::dirac(surf_sp) * L[surf_sp];
	Ms << NiHu::dirac(surf_sp) * M[surf_sp];
	Mts << NiHu::dirac(surf_sp) * Mt[surf_sp];
	Ds << NiHu::dirac(surf_sp) * D[surf_sp];
	Gxxs << NiHu::dirac(surf_sp) * Gxx[surf_sp];
	
	Lf << NiHu::dirac(field_sp) * L[surf_sp];
	Mf << NiHu::dirac(field_sp) * M[surf_sp];
	Mtf << NiHu::dirac(field_sp) * Mt[surf_sp];
	Df << NiHu::dirac(field_sp) * D[surf_sp];
}

