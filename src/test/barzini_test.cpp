#include "../bem/integral_operator.hpp"
#include "../library/poisson_kernel.hpp"

typedef tmp::vector<quad_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

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

	uMatrix elements(4, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		quad_1_elem::id, 4, 5, 8, 7;

	mesh_t msh(nodes, elements);

	auto f_sp = create_function_space_view(msh, field_option::isoparametric());

	dMatrix A(f_sp.get_num_dofs(), f_sp.get_num_dofs());
	A.setZero();
	
	auto b_op = create_integral_operator(poisson_G_kernel(), operator_option::non_local());
	eval_into(A, b_op(f_sp, f_sp, formalism::general()));

	std::cout << A << std::endl << std::endl;
	std::cout << b_op.get_kernel().get_num_evaluations() << std::endl;
	
	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / 4.0/M_PI;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.0)) << std::endl;

	return 0;
}

