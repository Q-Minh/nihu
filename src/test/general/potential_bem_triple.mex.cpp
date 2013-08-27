#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"
#include "library/poisson_singular_integrals.hpp"

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
	surf_elements <<
		tria_1_elem::id, 0, 1, 2;

	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _tria_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	dMatrix field_nodes(3,3);
	uMatrix field_elements(1,4);

	field_nodes <<
		0.0, 0.0, 1.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 1.0;
	field_elements <<
		tria_1_elem::id, 0, 1, 2;

	auto field_mesh = create_mesh(field_nodes, field_elements, _tria_1_tag());
	auto const &field_sp = constant_view(field_mesh);

	// field point system matrices

	auto n = surf_sp.get_num_dofs();
	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m,n); Lf.setZero();

	Lf << ( dirac(field_sp) * L[surf_sp] );

	std::cout << Lf << std::endl;
}

