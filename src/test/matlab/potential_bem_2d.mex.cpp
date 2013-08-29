#include "util/mex_matrix.hpp"
#include "bem/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 6 || nrhs < 4)
		throw("Too few input or output arguments");

	dMatrix surf_nodes(rhs[0]), surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _line_1_tag());
	auto const &trial_sp = constant_view(surf_mesh);
	auto const &test_sp = dirac(trial_sp);

	dMatrix field_nodes(rhs[2]), field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, _line_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	auto L = create_integral_operator(laplace_2d_SLP_kernel());
	auto M = create_integral_operator(laplace_2d_DLP_kernel());
	auto Mt = create_integral_operator(laplace_2d_DLPt_kernel());
	auto N = create_integral_operator(laplace_2d_HSP_kernel());
	auto I = identity_integral_operator();

	auto n = trial_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]), Mts(n, n, lhs[2]), Ns(n, n, lhs[3]);

	Ls << ( test_sp * L[trial_sp] );
	Ms << ( test_sp * M[trial_sp] ) + ( test_sp * (-.5*I)[trial_sp] );
	Mts << ( test_sp * Mt[trial_sp] ) + ( test_sp * (.5*I)[trial_sp] );
	Ns << ( test_sp * N[trial_sp] );

	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m, n, lhs[4]), Mf(m, n, lhs[5]);

	Lf << (field_sp * L[trial_sp]);
	Mf << (field_sp * M[trial_sp]);
}

