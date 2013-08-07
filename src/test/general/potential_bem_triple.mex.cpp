#include "bem/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

int main(void)
{
	// generating integral operators
	auto L = create_integral_operator(helmholtz_3d_SLP_kernel(1.0));

	// generating function spaces
	dMatrix surf_nodes(3,3);
	uMatrix surf_elements(1,4);

	surf_nodes <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0;
	surf_elements << tria_1_elem::id, 0, 1, 2;

	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _tria_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	// surface system matrices

	auto n = surf_sp.get_num_dofs();
	cMatrix Ls(n,n);

	Ls << ( dirac(surf_sp) * L[surf_sp] );

	std::cout << Ls << std::endl;
}

