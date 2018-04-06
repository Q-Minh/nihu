#include "library/helmholtz_kernel.hpp"
#include <gtest/gtest.h>

typedef Eigen::Matrix<double, 3, 1> x_t;

TEST(Helmholtz3dKernel, NormalDerivative)
{
	x_t x, y, nx, ny;
	x << 0., 0., 0.;
	y << 1., 0., 0.;
	nx << 3., 2., 1.;
	ny << 1., 2., 3.;
	nx = nx.normalized();
	ny = ny.normalized();
	
	// instantiate kernels
	NiHu::Helmholtz_3d_SLP_kernel G(k);
	NiHu::Helmholtz_3d_DLP_kernel Gy(k);
	NiHu::Helmholtz_3d_DLPt_kernel Gx(k);
	NiHu::Helmholtz_3d_SLP_kernel Gxy(k);
	
	// evaluate kernels
	std::complex<double> g = G(x, y);
	std::complex<double> gy = Gy(x, y, ny);
	std::complex<double> gx = Gx(x, y, nx);
	std::complex<double> gxy = Gxy(x, y, nx, ny);
	
	// approximate kernels
	double eps = 1e-2;
	std::complex<double> gy_num = (G(x, y+eps*ny) - G(x, y)) / eps;
	std::complex<double> gx_num = (G(x+eps*nx, y) - G(x, y)) / eps;
	std::complex<double> gxy_num = (G(x+eps*nx, y + eps*ny) - G(x+eps*nx, y)
		- G(x, y + eps*ny) + G(x, y)) / (eps*eps);
	
	// compute errors
	double err_y = std::abs(gy_num - gy) / std::abs(gy);
	double err_x = std::abs(gx_num - gy) / std::abs(gx);
	double err_xy = std::abs(gxy_num - gy) / std::abs(gxy);
	
	EXPECT_LE(err_y, 1e-5);
	EXPECT_LE(err_x, 1e-5);
	EXPECT_LE(err_xy, 1e-5);
}
