#include "core/weighted_residual.hpp"
#include "library/elastostatics_kernel.hpp"
#include "library/elastostatics_singular_integrals.hpp"
#include "library/lib_element.hpp"
#include "interface/read_off_mesh.hpp"

#include <fstream>
#include <sstream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

int main(int argc, char *argv[])
{
	auto mesh = read_off_mesh(argv[1], tria_1_tag(), quad_1_tag());
	auto const &v = constant_view(mesh, _3d());

	auto field = read_off_mesh(argv[2], quad_1_tag());
	auto const &w = constant_view(field, _3d());

	double nu = .33;
	auto U_op = create_integral_operator(elastostatics_3d_U_kernel(nu));
	auto T_op = create_integral_operator(elastostatics_3d_T_kernel(nu));

	dMatrix Us(v.get_num_dofs(), v.get_num_dofs());
	Us.setZero();
	Us << dirac(v) * U_op[v];

	dMatrix Ts(v.get_num_dofs(), v.get_num_dofs());
	Ts.setZero();
	Ts << dirac(v) * T_op[v];

	dMatrix Uf(w.get_num_dofs(), v.get_num_dofs());
	Uf.setZero();
	Uf << dirac(w) * U_op[v];

	dMatrix Tf(w.get_num_dofs(), v.get_num_dofs());
	Tf.setZero();
	Tf << dirac(w) * T_op[v];

	std::cout << "ready" << std::endl;

	std::ofstream o1;

	o1.open("Us.mtx", std::ofstream::out); o1 << Us; o1.close();
	o1.open("Uf.mtx", std::ofstream::out); o1 << Uf; o1.close();
	o1.open("Ts.mtx", std::ofstream::out); o1 << Ts; o1.close();
	o1.open("Tf.mtx", std::ofstream::out); o1 << Tf; o1.close();

	return 0;
}

