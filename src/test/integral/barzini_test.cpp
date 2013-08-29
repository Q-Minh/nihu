#include "bem/weighted_residual.hpp"
#include "library/poisson_kernel.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

template <class func_space>
void tester(func_space const &func_sp)
{
	int nDOF = func_sp.get_num_dofs();
	dMatrix A(nDOF, nDOF);
	A.setZero();

	auto L = create_integral_operator(poisson_3d_SLP_kernel());
	A << ( dirac(func_sp) * L[func_sp] );

	std::cout << "WR matrix:\n" << A << std::endl;
	std::cout << "sum of elements: " << A.sum() << std::endl;

	double anal = 32. * (std::log(1.+std::sqrt(2.))-(std::sqrt(2.)-1.)/3.) / 4./M_PI;
	std::cout << "log10 error = " << log10(std::abs(A.sum() / anal - 1.)) << std::endl;
}


auto build_mesh(void) ->
	decltype(create_mesh(dMatrix(), uMatrix(), _quad_1_tag(), _tria_1_tag()))
{
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

	uMatrix elements(5, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	return create_mesh(nodes, elements, _quad_1_tag(), _tria_1_tag());
}


int main(void)
{
	auto msh = build_mesh();

	std::cout << "Testing with constant field" << std::endl;
	tester(constant_view(msh));

	/*
	std::cout << "Testing with isoparametric field" << std::endl;
	tester(isoparametric_view(msh));
	*/

	return 0;
}

