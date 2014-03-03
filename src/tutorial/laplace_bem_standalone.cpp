#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix surf_nodes(6,3);
	surf_nodes <<
	0, 0, 0,
	1, 0, 0,
	2, 0, 0,
	0, 1, 0,
	1, 1, 0,
	2, 1, 0;

	uMatrix surf_elem(2,5);
	surf_elem <<
	quad_1_elem::id, 0, 1, 4, 3,
	quad_1_elem::id, 1, 2, 5, 4;

	auto surf_mesh = create_mesh(surf_nodes, surf_elem, quad_1_tag());

	auto const &surf_sp = constant_view(surf_mesh);

	int n = surf_sp.get_num_dofs();
	dMatrix L_surf(n, n), M_surf(n, n), Mt_surf(n, n), N_surf(n, n);
	L_surf.setZero();

	auto I = identity_integral_operator();
	auto L = create_integral_operator(laplace_3d_SLP_kernel());
	auto M = create_integral_operator(laplace_3d_DLP_kernel());
	auto Mt = create_integral_operator(laplace_3d_DLPt_kernel());
	auto N = create_integral_operator(laplace_3d_HSP_kernel());

	// conventional equations
	L_surf << dirac(surf_sp) * L[surf_sp];
	M_surf << dirac(surf_sp) * M[surf_sp]  +  dirac(surf_sp) * (-.5*I)[surf_sp];
	// hypersingular equations
	Mt_surf  << dirac(surf_sp) * Mt[surf_sp] +  dirac(surf_sp) * (.5*I)[surf_sp];
	N_surf  << dirac(surf_sp) * N[surf_sp];

	std::cout << "L:\t" << L_surf << std::endl;
	std::cout << "M:\t" << M_surf << std::endl;
}

