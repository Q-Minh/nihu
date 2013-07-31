#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

int main(void)
{
	// generating integral operators

	auto LM = create_integral_operator(
		create_couple_kernel(
			poisson_SLP_kernel(),
			poisson_DLP_kernel())
		);
	auto I = -.5 * identity_integral_operator();

	// generating function spaces
	dMatrix surf_nodes(4,3);
	uMatrix surf_elements(1,5);

	surf_nodes <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0,
		0.0, 1.0, 0.0;
	surf_elements << quad_1_elem::id, 0, 1, 2, 3;

	auto surf_mesh = create_mesh(surf_nodes, surf_elements, _quad_1_tag());
	auto const &surf_sp = constant_view(surf_mesh);

	// surface system matrices

	auto n = surf_sp.get_num_dofs();
	dMatrix Ls(n,n), Ms(n,n);

	(Ls, Ms) << ( surf_sp * LM[surf_sp] );
	Ms << ( surf_sp *  I[surf_sp] );

	// field point system matrices

	dMatrix field_nodes(4,3);
	uMatrix field_elements(1,5);

	field_nodes <<
		0.0, 0.0, 1.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 1.0,
		0.0, 1.0, 1.0;
	field_elements << quad_1_elem::id, 0, 1, 2, 3;


	auto field_mesh = create_mesh(field_nodes, field_elements, _quad_1_tag());
	auto const &field_sp = dirac(constant_view(field_mesh));

	auto m = field_sp.get_num_dofs();
	dMatrix Lf(m,n), Mf(m,n), Mf2(m,n);

	(Lf, Mf) << ( field_sp * LM[surf_sp] );

	std::cout << Mf;
}

