#include <boost/math/constants/constants.hpp>

#include "library/lib_element.hpp"
#include "core/mesh.hpp"
#include "library/elastostatics_kernel.hpp"
#include "library/elastostatics_singular_integrals.hpp"

#include "core/weighted_residual.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

/// \todo this mesh creating function should be extracted into a general module
NiHu::mesh<tmp::vector<NiHu::line_1_elem> >
CreateCircle(double R, unsigned N, float orient = 1.0)
{
	using namespace boost::math::double_constants;
	
	dMatrix nodes(N, 2);		// nodal locations
	for (unsigned i = 0; i < N; ++i)
	{
		double phi = orient * i * (2 * pi / N);
		nodes(i, 0) = R * cos(phi);
		nodes(i, 1) = R * sin(phi);
	}

	uMatrix elements(N, 3);		// element connections
	for (unsigned i = 0; i < N; ++i)
	{
		elements(i, 0) = NiHu::line_1_elem::id;
		elements(i, 1) = i;
		elements(i, 2) = (i + 1) % N;
	}

	return create_mesh(nodes, elements, NiHu::line_1_tag());
}

int main(void)
{
	float R = 1.0;
	auto	mesh = CreateCircle(R, 500);
	auto const &space = NiHu::constant_view(mesh, NiHu::_2d());
	auto const &tst_space = space;

	int n = space.get_num_dofs();
	std::cout << "number of DOFs: " << n << std::endl;
	dMatrix II(n, n);
	dMatrix LL(n, n);
	dMatrix MM(n, n);
	II.setZero();
	LL.setZero();
	MM.setZero();

	auto I = NiHu::identity_integral_operator();
	float nu = .3, mu = 1.0;
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
