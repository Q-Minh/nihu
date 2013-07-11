#include "../bem/integral_operator.hpp"
#include "../tmp/vector.hpp"
#include "../bem/mesh.hpp"
#include "../bem/function_space.hpp"
#include "../bem/projection.hpp"
#include "../bem/weighted_residual.hpp"

#include "../library/unit_kernel.hpp"
#include "../library/helmholtz_kernel.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
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

	auto trial = constant_view(msh);	// isoparametric
	auto test = constant_view(msh);
	int nDOF = test.get_num_dofs();
	dMatrix A(nDOF, nDOF);
	A.setZero();

	auto op = identity_integral_operator();
	auto wr = test * op[trial];
	wr.eval(A);

	std::cout << A << std::endl;

	return 0;
}

