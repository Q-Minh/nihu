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
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]) field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surface = create_mesh(surf_nodes, surf_elem, _quad_1_tag());
	auto field = create_mesh(field_nodes, field_elem, _quad_1_tag());
//! [Meshes]

//! [Function spaces]
	auto const &trial = isoparametric_view(surface);
	auto const &surf_test = dirac(constant_view(surface));
	auto const &field_test = dirac(constant_view(field));
//! [Function spaces]

//! [Kernel and weighted residual]
	double wave_number = *mxGetPr(rhs[4]);
	auto K = create_integral_operator(helmholtz_3d_SLP_kernel<double>(wave_number));

	cMatrix Z_trans(field_test.get_num_dofs(), trial.get_num_dofs(), lhs[0]);
	cMatrix Z_rad(surf_test.get_num_dofs(), trial.get_num_dofs(), lhs[1]);

	Z_trans << (field_test * K[trial]);
	Z_rad << (surf_test * K[trial]);
//! [Kernel and weighted residual]
}

