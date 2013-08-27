#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(9,3);
	nodes <<
		-1.0, -1.0, 0.0,
		 0.0, -1.0, 0.0,
		 1.0, -1.0, 0.0,

		-1.0,  0.0, 0.0,
		 0.0,  0.0, 0.0,
		 1.0,  0.0, 0.0,

		-1.0,  1.0, 0.0,
		 0.0,  1.0, 0.0,
		 1.0,  1.0, 0.0;

	uMatrix elements(5, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	auto msh = create_mesh(nodes, elements, _quad_1_tag(), _tria_1_tag());
	auto const &trial_sp = constant_view(msh);	// isoparametric
	auto const &test_sp = trial_sp;				// galerkin

	int nDOF = trial_sp.get_num_dofs();
	dMatrix A(nDOF, nDOF);
	A.setZero();

	auto L = create_integral_operator(poisson_3d_SLP_kernel());
	auto I = identity_integral_operator();

	A << (test_sp * L[trial_sp]) + (test_sp * (-.5*I)[trial_sp]);

	std::cout << A << std::endl << std::endl << A.sum() << std::endl;
	
	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / 4.0/M_PI;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.0)) << std::endl;

	return 0;
}

