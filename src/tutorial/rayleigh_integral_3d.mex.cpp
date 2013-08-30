//! [Header]
#include "bem/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "util/mex_matrix.hpp"

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;
//! [Header]

//! [Mex function]
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
//! [Mex function]
{
//! [Meshes]
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]), field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surface = create_mesh(surf_nodes, surf_elem, _quad_1_tag());
	auto field = create_mesh(field_nodes, field_elem, _quad_1_tag());
//! [Meshes]

//! [Function spaces]
	auto const &w = isoparametric_view(surface);
	auto const &dirac_s = dirac(constant_view(surface));
	auto const &dirac_f = dirac(constant_view(field));
//! [Function spaces]

//! [Kernel]
	double wave_number = *mxGetPr(rhs[4]);
	auto G = create_integral_operator(helmholtz_3d_SLP_kernel<double>(wave_number));
//! [Kernel]

//! [Matrices]
	cMatrix Z_rad(dirac_s.get_num_dofs(), w.get_num_dofs(), lhs[1]);
	cMatrix Z_trans(dirac_f.get_num_dofs(), w.get_num_dofs(), lhs[0]);
//! [Matrices]

//! [Weighted residual]
	Z_rad << (dirac_s * G[w]);
	Z_trans << (dirac_f * G[w]);
//! [Weighted residual]
}

