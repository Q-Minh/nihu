#include "../bem/integral_operator.hpp"
#include "../library/unit_kernel.hpp"

typedef tmp::vector<line_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(3,2);
	nodes <<
		-1.0, 0.0,
		 0.0, 0.0,
		 1.0, 0.0;

	uMatrix elements(2, 1+2);
	elements <<
		line_1_elem::id, 0, 1,
		line_1_elem::id, 1, 2;

	mesh_t msh(nodes, elements);

	auto f_sp = create_function_space_view(msh, field_option::constant());

	dMatrix A(f_sp.get_num_dofs(), f_sp.get_num_dofs());
	A.setZero();
	
	auto b_op = create_integral_operator(unit_kernel<space_2d>(), operator_option::local());
	A += b_op(f_sp, f_sp, formalism::full_dirac());

	std::cout << A << std::endl << std::endl;
	std::cout << b_op.get_kernel().get_num_evaluations() << std::endl;

	return 0;
}

