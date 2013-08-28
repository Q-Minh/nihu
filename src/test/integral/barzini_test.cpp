#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

template <class func_space>
void tester(function_space_base<func_space> const &func_sp)
{
	int nDOF = func_sp.get_num_dofs();
	dMatrix A(nDOF, nDOF);
	A.setZero();

	auto L = create_integral_operator(poisson_3d_SLP_kernel());
	A << ( func_sp.derived() * L[func_sp.derived()] );

	std::cout << A << std::endl << std::endl
		<< "sum of elements: " << A.sum() << std::endl;

	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / 4.0/M_PI;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.0)) << std::endl;
}


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

	std::cout << "Testing with constant field" << std::endl;
	tester(isoparametric_view(msh));

	std::cout << "Testing with isoparametric field" << std::endl;
	tester(constant_view(msh));

	return 0;
}

