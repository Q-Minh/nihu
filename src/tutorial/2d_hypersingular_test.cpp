#include "core/weighted_residual.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/lib_element.hpp"

int main(void)
{
	Eigen::Matrix<double, Eigen::Dynamic, 2> nodes(2, 2);
	nodes << -1.0, 0.0, 1.0, 0.0;
	
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> elements(1, 3);
	elements << NiHu::line_1_elem::id, 0, 1;
	
	auto mesh = NiHu::create_mesh(nodes, elements, NiHu::line_1_tag());
	auto const &space = NiHu::constant_view(mesh);
	size_t n = space.get_num_dofs();
	std::cout << "Number of DOF: " << n << std::endl;
	
	Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A(n, n);
	A.setZero();
	double k = 1.0;
	auto W = NiHu::create_integral_operator(NiHu::helmholtz_2d_HSP_kernel<double>(k));
	A << dirac(space) * W[space];
	
	return 0;
}
