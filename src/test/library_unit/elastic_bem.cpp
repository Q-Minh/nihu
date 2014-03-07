#include <iostream>
#include "core/weighted_residual.hpp"
#include "library/lib_element.hpp"
#include "library/elastostatics_kernel.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	dMatrix nodes(4,3);
	nodes <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0,
		0.0, 1.0, 0.0;
	uMatrix elements(1, 5);
	elements << quad_1_elem::id, 0, 1, 2, 3;
	auto mesh = create_mesh(nodes, elements, quad_1_tag());
	
	nodes <<
		0.0, 0.0, 2.0,
		1.0, 0.0, 2.0,
		1.0, 1.0, 2.0,
		0.0, 1.0, 2.0;
	auto field = create_mesh(nodes, elements, quad_1_tag());
	
	auto K = create_integral_operator(elastostatics_3d_U_kernel(.33));
	
	auto const &w = isoparametric_view(field, _3d());
	auto const &v = constant_view(mesh, _3d());
	
	dMatrix res(w.get_num_dofs(), v.get_num_dofs());

	res.setZero();
	res << dirac(w) * K[v];
	std::cout << res << std::endl;
	
	return 0;
}

