#include "nihu/library/helmholtz_kernel.hpp"
#include <gtest/gtest.h>

typedef Eigen::Matrix<double, 3, 1> x_t;

TEST(Helmholtz3dKernel, NormalDerivative)
{
	x_t x, y, nx, ny;
	x << 0., 0., 0.;
	y << 1., 1., 1.;
	nx << 1., 0., 1.;
	ny << 1., 1., 0.;
	nx = nx.normalized();
	ny = ny.normalized();
	
	// instantiate kernels
	double k = 3.5;
	NiHu::helmholtz_3d_SLP_kernel<double> G(k);
	NiHu::helmholtz_3d_DLP_kernel<double> Gy(k);
	NiHu::helmholtz_3d_DLPt_kernel<double> Gx(k);
	NiHu::helmholtz_3d_HSP_kernel<double> Gxy(k);
	
	// evaluate kernels
	std::complex<double> g = G(x, y);
	(void)g;
	std::complex<double> gy = Gy(x, y, ny);
	std::complex<double> gx = Gx(x, y, nx);
	std::complex<double> gxy = Gxy(x, y, nx, ny);
	
	// approximate kernels
	double eps = 1e-5;
	std::complex<double> gy_num = (G(x, y+eps*ny) - G(x, y)) / eps;
	std::complex<double> gx_num = (G(x+eps*nx, y) - G(x, y)) / eps;
	std::complex<double> gxy_num = (G(x+eps*nx, x_t(y + eps*ny)) - G(x+eps*nx, y)
		- G(x, y + eps*ny) + G(x, y)) / (eps*eps);
	
	// compute errors
	double err_y = std::abs(gy_num - gy) / std::abs(gy);
	double err_x = std::abs(gx_num - gx) / std::abs(gx);
	double err_xy = std::abs(gxy_num - gxy) / std::abs(gxy);

	std::cout << err_y << '\n' << err_x << '\n' << err_xy << '\n';
	
	EXPECT_LE(err_y, 1e-4);
	EXPECT_LE(err_x, 1e-4);
	EXPECT_LE(err_xy, 1e-4);
}
