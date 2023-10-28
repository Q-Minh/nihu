#include <mex.h>

#include "nihu/library/quad_1_gauss_field.hpp"
#include "nihu/core/weighted_residual.hpp"
#include "nihu/util/mex_matrix.hpp"
#include "nihu/library/helmholtz_kernel.hpp"
#include "nihu/library/helmholtz_nearly_singular_integrals.hpp"
#include "nihu/library/helmholtz_singular_integrals.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef NiHu::mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	// input arguments
	if (nrhs != 5)
		mexErrMsgTxt("Usage: mex(surf_nodes, surf_elem, field_nodes, field_elem, k)");
	dMatrix
 		surf_nodes(rhs[0]), surf_elem(rhs[1]),	
		field_nodes(rhs[2]), field_elem(rhs[3]);
	double k = *mxGetPr(rhs[4]);
	
	// assemble field matrix
	size_t nElements = surf_elem.rows();
	uMatrix fields(nElements, 1+4+4);
	for (size_t e = 0; e < nElements; ++e)
	{
		fields(e,0) = NiHu::quad_1_gauss_field::id;
		for (size_t c = 0; c < 4; ++c)
			fields(e,c+1) = surf_elem(e,c+1);
		for (size_t c = 0; c < 4; ++c)
			fields(e,c+1+4) = unsigned(4*e+c);
	}
	
	// create function space
	auto surf_sp = NiHu::create_function_space(surf_nodes, fields, NiHu::quad_1_gauss_field_tag());
		
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem,
		NiHu::tria_1_tag(), NiHu::quad_1_tag());
	
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field_mesh));

	size_t n = surf_sp.get_num_dofs();
	size_t m = field_sp.get_num_dofs();
	cMatrix
		L_surf(n, n, lhs[0]), M_surf(n, n, lhs[1]),
		Mt_surf(n, n, lhs[2]), D_surf(n, n, lhs[3]),
		L_field(m, n, lhs[4]), M_field(m, n, lhs[5]);
	
	auto L = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(k));
	auto M = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));
	auto Mt = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLPt_kernel<double>(k));
	auto D = NiHu::create_integral_operator(NiHu::helmholtz_3d_HSP_kernel<double>(k));
	
// #pragma omp parallel sections
// {
	// #pragma omp section
	// {
		L_surf << NiHu::dirac(surf_sp) * L[surf_sp];
	// }
	// #pragma omp section
	// {
	M_surf << NiHu::dirac(surf_sp) * M[surf_sp];
	// }
	// #pragma omp section
	// {
	Mt_surf << NiHu::dirac(surf_sp) * Mt[surf_sp];
	// }
	// #pragma omp section
	// {
	D_surf << NiHu::dirac(surf_sp) * D[surf_sp];
	// }

	// #pragma omp section
	// {
	L_field << field_sp * L[surf_sp];
	// }
	// #pragma omp section
	// {
	M_field << field_sp * M[surf_sp];
	// }
// }

}

