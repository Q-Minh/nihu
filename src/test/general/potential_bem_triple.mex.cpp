#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

int main(void)
{
	// generating integral operators

	auto L = create_integral_operator(poisson_3d_SLP_kernel());

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
	dMatrix Ls(n,n);

	Ls << ( dirac(surf_sp) * L[surf_sp] );

	std::cout << Ls << std::endl;
}

