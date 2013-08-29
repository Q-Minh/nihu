//! [Includes]
#include "bem/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
//! [Includes]

//! [Typedefs]
// dynamically resizable double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// dynamically resizable unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
//! [Typedefs]

//! [Test]
template <class func_space>
void tester(func_space const &w)
{
	// compute number of DOF and allocate result matrix
	int nDOF = w.get_num_dofs();
	dMatrix A(nDOF, nDOF);
	A.setZero();

	// create integral operator from kernel and perform weighted double integral
	auto K = create_integral_operator(laplace_3d_SLP_kernel());
	A << ( w * K[w] );

	// Display matrix elements and their sum
	std::cout << "WR matrix:\n" << A << std::endl;
	std::cout << "sum of elements: " << A.sum() << std::endl;

	// Compare to analytical solution
	double anal = 32. * (std::log(1.+std::sqrt(2.))-(std::sqrt(2.)-1.)/3.) / 4./M_PI;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.)) << std::endl;
}
//! [Test]


//! [Mesh generation]
auto build_mesh(void) ->
	decltype(create_mesh(dMatrix(), uMatrix(), _tria_1_tag(), _quad_1_tag()))
{
	// nodal coordinates in 3D, 9 nodes
	dMatrix nodes(9,3);
	nodes <<
		-1., -1., 0.,
		 0., -1., 0.,
		 1., -1., 0.,
		//--------------
		-1.,  0., 0.,
		 0.,  0., 0.,
		 1.,  0., 0.,
		//--------------
		-1.,  1., 0.,
		 0.,  1., 0.,
		 1.,  1., 0.;

	// element nodal indices
	uMatrix elements(5, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	// create a mesh consising of tria and quad elements
	return create_mesh(nodes, elements, _tria_1_tag(), _quad_1_tag());
}
//! [Mesh generation]


//! [Main]
int main(void)
{
	// call the mesh generating function
	auto mesh = build_mesh();

	// create a pieceweise constant function space and call the tester
	std::cout << "Testing with constant field" << std::endl;
	tester(constant_view(mesh));

	// create a pieceweise linear function space and call the tester
	std::cout << "Testing with isoparametric field" << std::endl;
	tester(isoparametric_view(mesh));

	return 0;
}
//! [Main]

