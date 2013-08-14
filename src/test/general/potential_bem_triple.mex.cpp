#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"
#include "library/poisson_singular_integrals.hpp"

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

int main(void)
{
	// generating integral operators
	auto L = create_integral_operator(poisson_2d_SLP_kernel());
	auto M = create_integral_operator(poisson_2d_DLP_kernel());

	// generating function spaces
	dMatrix surf_nodes(3,2);
	uMatrix surf_elements(2,3);

	surf_nodes <<
		0.0, 0.0,
		1.0, 0.0,
		2.0, 0.0;
	surf_elements <<
		line_1_elem::id, 0, 1,
		line_1_elem::id, 1, 2;

	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _line_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	// surface system matrices

	auto n = surf_sp.get_num_dofs();
	dMatrix Ls(n,n); Ls.setZero();
	dMatrix Ms(n,n); Ms.setZero();

	Ls << ( dirac(surf_sp) * L[surf_sp] );
	Ms << ( dirac(surf_sp) * M[surf_sp] );

	std::cout << Ls << std::endl;
	std::cout << Ms << std::endl;
}

