#include "core/weighted_residual.hpp"
#include "library/elastostatics_kernel.hpp"
#include "library/lib_element.hpp"
#include "interface/read_off_mesh.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

int main(int argc, char *argv[])
{
	auto mesh = read_off_mesh(argv[1], tria_1_tag(), quad_1_tag());
	auto const &v = constant_view(mesh, _3d());

	auto field = read_off_mesh(argv[2], tria_1_tag(), quad_1_tag());
	auto const &w = constant_view(field, _3d());

	auto L_op = create_integral_operator(elastostatics_3d_U_kernel(.33));

	dMatrix L(w.get_num_dofs(), v.get_num_dofs());
	L.setZero();
	L << dirac(w) * L_op[v];

	std::cout << L << std::endl;

	return 0;
}

