#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"
#include "library/lib_element.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
typedef NiHu::mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]), field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::line_1_tag());
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem, NiHu::quad_1_volume_tag());

	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field_mesh));

	int n = surf_sp.get_num_dofs();
	int m = field_sp.get_num_dofs();
	cMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]), Ws(n, n, lhs[2]), Lf(m, n, lhs[3]), Mf(m, n, lhs[4]);

	double k = *mxGetPr(rhs[4]);
	auto L = NiHu::create_integral_operator(NiHu::helmholtz_2d_SLP_kernel<double>(k));
	auto M = NiHu::create_integral_operator(NiHu::helmholtz_2d_DLP_kernel<double>(k));
	auto W = NiHu::create_integral_operator(NiHu::helmholtz_2d_HSP_kernel<double>(k));

	Ls << dirac(surf_sp) * L[surf_sp]; 
	Ms << dirac(surf_sp) * M[surf_sp];
	Ws << dirac(surf_sp) * W[surf_sp];
	Lf  << field_sp * L[surf_sp];
	Mf  << field_sp * M[surf_sp];
}

