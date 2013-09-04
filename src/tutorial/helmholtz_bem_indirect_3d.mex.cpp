//![Header]
#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef mex::real_matrix<double> dMatrix;
//![Header]

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
//![Meshes]
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elem, _tria_1_tag());
	dMatrix field_nodes(rhs[2]), field_elem(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elem, _quad_1_tag());
//![Meshes]
	
//! [Function spaces]
	auto const &surf_sp = constant_view(surf_mesh);
	auto const &field_sp = dirac(constant_view(field_mesh));
//! [Function spaces]

//! [Matrices]
	int n = surf_sp.get_num_dofs();
	int m = field_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]), Mts(n, n, lhs[1]), Ns(n, n, lhs[2]),
		Lf(m, n, lhs[3]), Mf(m, n, lhs[4]);
//! [Matrices]
	
//! [Integral operators]
	auto I = identity_integral_operator();
	auto L = create_integral_operator(laplace_3d_SLP_kernel());
	auto M = create_integral_operator(laplace_3d_DLP_kernel());
	auto Mt = create_integral_operator(laplace_3d_DLPt_kernel());
	auto N = create_integral_operator(laplace_3d_HSP_kernel());
//! [Integral operators]

//! [System matrices]
	Ls << dirac(surf_sp) * L[surf_sp];
	Mts << dirac(surf_sp) * Mt[surf_sp] + dirac(surf_sp) * (.5*I)[surf_sp];
	Ns << dirac(surf_sp) * N[surf_sp];
	Lf << field_sp * L[surf_sp];
	Mf << field_sp * M[surf_sp];
//! [System matrices]
}

