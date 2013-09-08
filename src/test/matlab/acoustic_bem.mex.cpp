#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/helmholtz_kernel.hpp"

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elem, _quad_1_tag(), _tria_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);
	int n = surf_sp.get_num_dofs();

	dMatrix field_nodes(rhs[2]), field_elem(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elem, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));
	int m = field_sp.get_num_dofs();

	cMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]), Lf(m, n, lhs[2]), Mf(m, n, lhs[3]);
	
	double k = *mxGetPr(rhs[4]);
	auto I = identity_integral_operator();
	auto Op = create_integral_operator(
		helmholtz_3d_SLP_kernel<double>(k),
		helmholtz_3d_DLP_kernel<double>(k)
	);

	(Ls, Ms) << dirac(surf_sp) * Op[surf_sp];
	Ms << dirac(surf_sp) * (-.5*I)[surf_sp];
	(Lf, Mf)  << field_sp * Op[surf_sp];
}

