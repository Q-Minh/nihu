#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_singular_integrals.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(4,2);
	uMatrix elements(4,3);

	nodes <<
		0.0, 0.0,
		1.0, 0.0,
		1.0, 1.0,
		0.0, 1.0;
	elements <<
		line_1_elem::id, 0, 1,
		line_1_elem::id, 1, 2,
		line_1_elem::id, 2, 3,
		line_1_elem::id, 3, 0;

	auto mesh = create_mesh(nodes, elements, line_1_tag());
	auto const &fspace = constant_view(mesh);
	auto nDof = fspace.get_num_dofs();

	dMatrix L(nDof, nDof), M(nDof, nDof);
	L.setZero(); M.setZero();

	auto L_op = create_integral_operator(laplace_2d_SLP_kernel());
	auto M_op = create_integral_operator(laplace_2d_DLP_kernel());

	L << fspace * L_op[fspace];
	M << fspace * M_op[fspace];

	std::cout << L << std::endl;
	std::cout << M << std::endl;

	return 0;
}
