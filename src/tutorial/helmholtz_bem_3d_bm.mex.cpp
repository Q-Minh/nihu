#include <mex.h>

#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
typedef NiHu::mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::tria_1_tag(), NiHu::quad_1_tag());
	
	auto const &surf_sp = NiHu::constant_view(surf_mesh);

	size_t n = surf_sp.get_num_dofs();
	cMatrix L_surf(n, n, lhs[0]), M_surf(n, n, lhs[1]),
		Mt_surf(n, n, lhs[2]), N_surf(n, n, lhs[3]);
	
	double k = *mxGetPr(rhs[2]);
	auto L = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(k));;
	auto M = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));
	auto Mt = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLPt_kernel<double>(k));
	auto N = NiHu::create_integral_operator(NiHu::helmholtz_3d_HSP_kernel<double>(k));

	L_surf << dirac(surf_sp) * L[surf_sp];
	M_surf << dirac(surf_sp) * M[surf_sp];
	Mt_surf << dirac(surf_sp) * Mt[surf_sp];
	N_surf << dirac(surf_sp) * N[surf_sp];
}

