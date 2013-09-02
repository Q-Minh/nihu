//![Header]
#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;
//![Header]

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
//![Mesh]
	dMatrix
 		surf_nodes(rhs[0]), surf_elem(rhs[1]),	
		chief_nodes(rhs[2]), chief_elem(rhs[3]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elem, _tria_1_tag());
	auto chief_mesh = create_mesh(chief_nodes, chief_elem, _quad_1_tag());
//![Mesh]
	
//! [Function spaces]
	auto const &surf_sp = constant_view(surf_mesh);
	auto const &chief_sp = constant_view(chief_mesh);
//! [Function spaces]

//! [Matrices]
	int n = surf_sp.get_num_dofs();
	int m = chief_sp.get_num_dofs();
	cMatrix
		Ls(n, n, lhs[0]), Ms(n, n, lhs[1]),
		Ls_chief(m, n, lhs[2]), Ms_chief(m, n, lhs[3]),
		Mts(n, n, lhs[4]), Ns(n, n, lhs[5]);
//! [Matrices]
	
//! [Integral operators]
	double k = *mxGetPr(rhs[4]);
	auto I = identity_integral_operator();
	auto L = create_integral_operator(helmholtz_3d_SLP_kernel<double>(k));
	auto M = create_integral_operator(helmholtz_3d_DLP_kernel<double>(k));
	auto Mt = create_integral_operator(helmholtz_3d_DLPt_kernel<double>(k));
	auto N = create_integral_operator(helmholtz_3d_HSP_kernel<double>(k));
//! [Integral operators]

//! [System matrices]
	// plain
	Ls << dirac(surf_sp) * L[surf_sp];
	Ms << dirac(surf_sp) * M[surf_sp]  +  dirac(surf_sp) * (-.5*I)[surf_sp];
	// CHIEF
	Ls_chief << dirac(chief_sp) * L[surf_sp];
	Ms_chief << dirac(chief_sp) * M[surf_sp];
	// Burton-Miller
	Mts  << dirac(surf_sp) * Mt[surf_sp] +  dirac(surf_sp) * (.5*I)[surf_sp];
	Ns  << dirac(surf_sp) * N[surf_sp];
//! [System matrices]
}

