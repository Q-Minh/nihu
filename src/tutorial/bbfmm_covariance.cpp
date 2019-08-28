#include "library/covariance_kernel.hpp"
#include "fmm/black_box_fmm.hpp"

int main()
{
	typedef NiHu::covariance_kernel<NiHu::space_3d<> > kernel_t;

	double sigma = 1.;
	double length = 1.;
	kernel_t kernel(sigma, length);

	typedef NiHu::fmm::black_box_fmm<kernel_t> bbfmm_t;
	bbfmm_t bbfmm(kernel);

	auto p2p = bbfmm.create_p2p();
	auto p2m = bbfmm.create_p2m();
	auto m2l = bbfmm.create_m2l();
	auto l2p = bbfmm.create_l2p();

	typedef kernel_t::location_t location_t;
	location_t x, y;
	x << 0, 0, 0;
	y << 10, 0, 0;

	location_t X, Y;
	X << .5, 0, 0;
	Y << 9.5, 0, 0;

	typedef bbfmm_t::cluster_t cluster_t;
	cluster_t CX, CY;
	CX.set_bounding_box(cluster_t::bounding_box_t(X, 2.0));
	CY.set_bounding_box(cluster_t::bounding_box_t(Y, 2.0));
	CX.set_chebyshev_order(5);
	CY.set_chebyshev_order(5);

	auto decomp = (l2p(x, CX) * m2l(CX, CY) * p2m(CY, y)).eval();
	auto anal = p2p(x, y);

	std::cout << anal << std::endl;
	std::cout << decomp << std::endl;

	return 0;
}
