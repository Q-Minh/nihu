/**
 * @file laplace_2d_coll_mex.mex.cpp
 * @brief Conventional 2D BEM for the Laplace equation using collocational formalism
 * @details
 * @ingroup app_laplace
 */

#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/lib_element.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]), field_nodes(rhs[2]), field_elem(rhs[3]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::line_1_tag(), NiHu::line_2_tag());
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem, NiHu::tria_1_volume_tag(), NiHu::quad_1_volume_tag());

	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field_mesh));

	size_t n = surf_sp.get_num_dofs();
	size_t m = field_sp.get_num_dofs();
	dMatrix Ls(n, n, lhs[0]), Ms(n, n, lhs[1]), Lf(m, n, lhs[2]), Mf(m, n, lhs[3]);

	auto I = NiHu::identity_integral_operator();
	auto L = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto M = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());

	Ls << NiHu::dirac(surf_sp) * L[surf_sp]; 
	Ms << NiHu::dirac(surf_sp) * M[surf_sp] 
		+ NiHu::dirac(surf_sp) * (-.5*I)[surf_sp];
	Lf  << field_sp * L[surf_sp];
	Mf  << field_sp * M[surf_sp];
}

