#include "../bem/integral_operator.hpp"
#include "../bem/weighted_residual.hpp"
#include "../library/poisson_kernel.hpp"

typedef tmp::vector<quad_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

typedef function_space_view<mesh_t, isoparametric_field> function_space_t;

// double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	// define nodal coordinates of our mesh
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

	// define elements of our mesh
	uMatrix elements(4, 1+4);
	elements <<
	//  field id         nodes
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		quad_1_elem::id, 4, 5, 8, 7;

	// instantiate the mesh and the function space
	mesh_t msh(nodes, elements);
	function_space_t fsp(msh);

	// allocate and clear the result matrix
	dMatrix A(fsp.get_num_dofs(), fsp.get_num_dofs());
	A.setZero();
	
	// evaluate the operator into the matrix
	auto bound_op = boundary_operator(poisson_G_kernel(), non_local());
	eval_into(A, bound_op(fsp, fsp, formalism::general()) );

	// print the result matrix
	std::cout << A << std::endl << std::endl;
	// and the number of kernel evaluations
	std::cout << bound_op.get_kernel().get_num_evaluations() << std::endl;
	
	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / (4.0 * M_PI);
	std::cout << "numer: " << A.sum() << std::endl;
	std::cout << "anal: " << anal << std::endl;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.0)) << std::endl;

	return 0;
}

