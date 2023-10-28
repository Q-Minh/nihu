#include <boost/math/constants/constants.hpp>

#include "nihu/core/mesh.hpp"
#include "nihu/library/elastostatics_kernel.hpp"
#include "nihu/library/elastostatics_singular_integrals.hpp"
#include "nihu/library/lib_element.hpp"
#include "nihu/library/lib_mesh.hpp"

#include "nihu/core/weighted_residual.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	float R = 1.0;
	size_t N = 500;
	auto mesh = NiHu::create_line_1_circle_mesh(R, N);
	auto const &space = NiHu::constant_view(mesh, NiHu::field_dimension::_2d());
	auto const &tst_space = space;

	size_t n = space.get_num_dofs();
	std::cout << "number of DOFs: " << n << std::endl;
	dMatrix II(n, n);
	dMatrix LL(n, n);
	dMatrix MM(n, n);
	II.setZero();
	LL.setZero();
	MM.setZero();

	auto I = NiHu::identity_integral_operator();
	float nu = .3f, mu = 1.0f;
	auto L = NiHu::create_integral_operator(NiHu::elastostatics_2d_U_kernel(nu, mu));
	auto M = NiHu::create_integral_operator(NiHu::elastostatics_2d_T_kernel(nu, mu));

	II << tst_space * I[space];
	LL << tst_space * L[space];
	MM << tst_space * M[space];

	std::cout << II.topLeftCorner(10, 10) << std::endl << std::endl;
	std::cout << LL.topLeftCorner(10, 10) << std::endl << std::endl;
	std::cout << MM.topLeftCorner(10, 10) << std::endl;

	return 0;
}
