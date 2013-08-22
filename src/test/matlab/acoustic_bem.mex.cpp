#include "util/mex_matrix.hpp"
#include "bem/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"

#include <chrono>

typedef mex::real_matrix<double> dMatrix;
typedef mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs < 10 || nrhs < 5)
		throw("Too few input or output arguments");

	// create integral operators

	double k = *mxGetPr(rhs[4]);
	auto L = create_integral_operator(helmholtz_3d_SLP_kernel<double>(k));
	auto M = create_integral_operator(helmholtz_3d_DLP_kernel<double>(k));
	auto Mt = create_integral_operator(helmholtz_3d_DLPt_kernel<double>(k));
	auto N = create_integral_operator(helmholtz_3d_HSP_kernel<double>(k));
	auto I = identity_integral_operator();

	// create surface function spaces

	dMatrix surf_nodes(rhs[0]), surf_elements(rhs[1]);
	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag(), _tria_1_tag());
	auto const &trial_sp = constant_view(surf_mesh);
	auto const &test_sp = dirac(trial_sp);

	// generate field point system matrices

	auto n = trial_sp.get_num_dofs();
	cMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]), Mts(n, n, lhs[2]), Ns(n, n, lhs[3]);

	Ls << ( test_sp * L[trial_sp] );
	Ms << ( test_sp * M[trial_sp] ) + ( test_sp * (-.5*I)[trial_sp] );
	Mts << ( test_sp * Mt[trial_sp] ) + ( test_sp * (.5*I)[trial_sp] );
	Ns << ( test_sp * N[trial_sp] );

	// create field point function spaces

	dMatrix field_nodes(rhs[2]), field_elements(rhs[3]);
	auto field_mesh = create_mesh(field_nodes, field_elements, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	// generate field point system matrices

	auto m = field_sp.get_num_dofs();
	cMatrix Lf(m, n, lhs[4]), Mf(m, n, lhs[5]), Mtf(m, n, lhs[6]), Nf(m, n, lhs[7]);

	// evaluate radiation matrices with separate kernel evaluation

	auto start = std::chrono::steady_clock::now();	// start timer
	Lf  << ( field_sp * L [trial_sp] );
	Mf  << ( field_sp * M [trial_sp] );
	Mtf << ( field_sp * Mt[trial_sp] );
	Nf  << ( field_sp * N [trial_sp] );
	auto stop = std::chrono::steady_clock::now();	// stop timer
	dMatrix dur_separate(1, 1, lhs[8]);
	dur_separate(0,0) = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

	// evaluate radiation matrices again with couple kernel evaluation

	auto CoupleOp = create_integral_operator(
		helmholtz_3d_SLP_kernel<double>(k),
		helmholtz_3d_DLP_kernel<double>(k),
		helmholtz_3d_DLPt_kernel<double>(k),
		helmholtz_3d_HSP_kernel<double>(k)
	);
	start = std::chrono::steady_clock::now();	// start timer
	create_couple(Lf, Mf, Mtf, Nf)  << ( field_sp * CoupleOp [trial_sp] );
	stop = std::chrono::steady_clock::now();	// stop timer
	dMatrix dur_couple(1, 1, lhs[9]);
	dur_couple(0,0) = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
}

