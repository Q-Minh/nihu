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
//! [Surface mesh]
	auto surf_nodes = dMatrix(rhs[0]);
	auto surf_elem = dMatrix(rhs[1]);
	auto surface = create_mesh(surf_nodes, surf_elem, _quad_1_tag());
//! [Surface mesh]

//! [Field point mesh]
	auto field_nodes = dMatrix(rhs[2]);
	auto field_elem = dMatrix(rhs[3]);
	auto field = create_mesh(field_nodes, field_elem, _quad_1_tag());
//! [Field point mesh]
	
//! [Function spaces]
	auto const &trial = isoparametric_view(surface);
	auto const &test = dirac(constant_view(field));
//! [Function spaces]

//! [Kernel and weighted residual]
	double wave_number = *mxGetPr(rhs[4]);
	auto K = create_integral_operator(helmholtz_3d_SLP_kernel<double>(wave_number));
	cMatrix Z(test.get_num_dofs(), trial.get_num_dofs(), lhs[0]);
	Z << (test * K[trial]);
//! [Kernel and weighted residual]
}

