#include "core/weighted_residual.hpp"
#include "library/elastostatics_kernel.hpp"
#include "library/lib_element.hpp"
#include "interface/read_off_mesh.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

int main(int argc, char *argv[])
{
	auto mesh = read_off_mesh(argv[1], tria_1_tag());
	auto const &v = constant_view(mesh, _3d());

	auto field = read_off_mesh(argv[2], quad_1_tag());
	auto const &w = constant_view(field, _3d());

	double nu = .33;
	auto U_op = create_integral_operator(elastostatics_3d_U_kernel(nu));
	auto T_op = create_integral_operator(elastostatics_3d_T_kernel(nu));

	dMatrix U(w.get_num_dofs(), v.get_num_dofs());
	U.setZero();
	dMatrix T(w.get_num_dofs(), v.get_num_dofs());
	T.setZero();
	U << dirac(w) * U_op[v];
	T << dirac(w) * T_op[v];

	std::cout << U << std::endl;
	std::cout << "\n" << std::endl;
	std::cout << T << std::endl;

	return 0;
}

