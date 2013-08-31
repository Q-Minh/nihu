//! [Includes]
#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
//! [Includes]

//! [Typedefs]
// dynamically resizeable double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// dynamically resizeable unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
//! [Typedefs]

//! [Test]
template <class func_space>
void tester(func_space const &w)
{
	// compute number of DOF and allocate result matrix
	int nDOF = w.get_num_dofs();
	dMatrix I(nDOF, nDOF);
	I.setZero();

	// create integral operator from kernel and perform weighted double integral
	auto K = create_integral_operator(laplace_3d_SLP_kernel());
	I <<  w * K[w];

	// Display matrix elements and their sum
	std::cout << "WR matrix:\n" << I << std::endl;
	std::cout << "sum of elements: " << I.sum() << std::endl;

	// Compare to analytical solution
	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / 4.0/M_PI;
	std::cout << "log10 error = " << std::log10(std::abs(I.sum() / anal - 1.0)) << std::endl;
}
//! [Test]


//! [Main]
int main(void)
{
	// nodal coordinates in 3D, 9 nodes
	dMatrix nodes(4,3);
	nodes <<
		-1.0, -1.0, 0.0,
		 1.0, -1.0, 0.0,
		-1.0,  1.0, 0.0,
		 1.0,  1.0, 0.0;

	// element nodal indices
	uMatrix elements(1, 1+4);
	elements << quad_1_elem::id, 0, 1, 3, 2;

	// create the mesh
	auto msh = create_mesh(nodes, elements, _quad_1_tag());

	// create a piecewise constant function space and call the tester
	std::cout << "Testing with constant field" << std::endl;
	tester(constant_view(msh));

	// create a piecewise linear function space and call the tester
	std::cout << "Testing with isoparametric field" << std::endl;
	tester(isoparametric_view(msh));

	return 0;
}
//! [Main]

