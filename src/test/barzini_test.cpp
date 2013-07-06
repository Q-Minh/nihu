#include "../bem/integral_operator.hpp"
#include "../library/poisson_kernel.hpp"
#include "../library/unit_kernel.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector_t;
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

	uMatrix elements(5, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	mesh_t msh(nodes, elements);

	auto trial_sp = isoparametric_view(msh);	// isoparametric
	auto test_sp = trial_sp;				// collocational

	int nDOF = trial_sp.get_num_dofs();
	dMatrix A(nDOF, nDOF);
	A.setZero();

	auto b_op = nonlocal(poisson_G_kernel());
	auto id_op = local(unit_kernel<space_3d>());

	A += id_op(test_sp, trial_sp);

	std::cout << A << std::endl << std::endl << A.sum() << std::endl;
	std::cout << b_op.get_kernel().get_num_evaluations() << std::endl;
	
	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / 4.0/M_PI;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.0)) << std::endl;

	return 0;
}

