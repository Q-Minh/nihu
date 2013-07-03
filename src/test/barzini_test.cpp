#include "../bem/weighted_residual.hpp"
#include "../library/poisson_kernel.hpp"

typedef tmp::vector<quad_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;
typedef function_space_view<mesh_t, isoparametric_field> function_space_t;
typedef weighted_residual<GENERAL_FORMULATION, poisson_G_kernel, function_space_t> wr_t;

// double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	// instantiate kernel
	poisson_G_kernel kernel;

	// we define the nodal coordinates of our mesh
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

	// we define the fields of our function space
	uMatrix elements(4, 1+4);
	elements <<
	//  field id         nodes
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		quad_1_elem::id, 4, 5, 8, 7;

	// we instantiate the mesh using nodes and fields
	mesh_t msh(nodes, elements);
	// and the function space using our mesh
	function_space_t fsp(msh);

	// we allocate and clear the result matrix
	dMatrix A(fsp.get_num_dofs(), fsp.get_num_dofs());
	A.setZero();
	
	// and evaluate the weighed residual into our result matrix
	wr_t::eval(A, kernel, fsp, fsp);

	// we print the result matrix
	std::cout << A << std::endl << std::endl;
	// and the number of kernel evaluations
	std::cout << kernel.get_num_evaluations() << std::endl;
	
	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / 4.0/M_PI;
	std::cout << "numer: " << A.sum() << std::endl;
	std::cout << "anal: " << anal << std::endl;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.0)) << std::endl;

	return 0;
}

