#include "bem/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;


int main(void)
{
//! [Mesh]
	dMatrix nodes(4,3);
	nodes << 0., 0., 0., /*|*/ 1., 0., 0., /*|*/ 1., 1., 0., /*|*/ 0., 1., 0.;
	uMatrix elements(1, 1+4);
	elements <<	quad_1_elem::id, 0, 1, 2, 3;
	auto msh = create_mesh(nodes, elements, _quad_1_tag());
//! [Mesh]

//! [Function spaces]
	auto const &trial = constant_view(msh);
	auto const &test = dirac(trial);
//! [Function spaces]

//! [Weighted residual]
	dMatrix A(1, 1);
	A.setZero();

	auto K = create_integral_operator(laplace_3d_SLP_kernel());
	A << ( test * K[trial] );
//! [Weighted residual]

//! [Results]
	std::cout << "WR matrix: " << A << std::endl;

	double anal = std::log(1.+std::sqrt(2.)) / M_PI;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.)) << std::endl;
//! [Results]

	return 0;
}

