/**
 * @file laplace_2d_coll_mex.mex.cpp
 * @brief Conventional 2D BEM for the Laplace equation using collocational formalism
 * @details
 * @ingroup app_laplace
 */

#define NIHU_DEBUGGING 1

#include "core/weighted_residual.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/lib_element.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

int main(int argc, char *argv[])
{
	std::string surf_name(argv[1]);
	std::string field_name(argv[2]);

	auto surf_mesh = NiHu::read_off_mesh(surf_name, NiHu::line_2_tag());
	auto field_mesh = NiHu::read_off_mesh(field_name, NiHu::quad_1_volume_tag());

	auto const &surf_sp = NiHu::constant_view(surf_mesh);
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field_mesh));

	size_t n = surf_sp.get_num_dofs();
	size_t m = field_sp.get_num_dofs();
	dMatrix Ls = dMatrix::Zero(n, n);
	dMatrix Ms = dMatrix::Zero(n, n);
	dMatrix Lf = dMatrix::Zero(m, n);
	dMatrix Mf = dMatrix::Zero(m, n);

	auto I = NiHu::identity_integral_operator();
	auto L = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());
	auto M = NiHu::create_integral_operator(NiHu::laplace_2d_DLP_kernel());

	Ms << NiHu::dirac(surf_sp) * M[surf_sp]
		+ NiHu::dirac(surf_sp) * (-.5 * I)[surf_sp];
	Ls << NiHu::dirac(surf_sp) * L[surf_sp];
	Lf  << field_sp * L[surf_sp];
	Mf  << field_sp * M[surf_sp];

	return 0; 
}

