#include <boost/math/constants/constants.hpp>

#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/lib_element.hpp"
#include "nihu/core/inverse_mapping.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iomanip>

template <class elem_t, class nset_t>
Eigen::Matrix<double, 4, 1>
stokes_hsp_integral(
	elem_t const &elem,
	typename elem_t::x_t const &x,
	typename elem_t::x_t const &nx,
	typename elem_t::x_t const &y0,
	typename nset_t::shape_t const &N0,
	Eigen::Matrix<double, 4, 3> const &gradN0
)
{
	using namespace boost::math::double_constants;

	typedef typename elem_t::x_t x_t;
	typedef typename elem_t::xi_t xi_t;
	typedef typename elem_t::lset_t lset_t;
	typedef typename elem_t::domain_t domain_t;

	Eigen::Matrix<double, 4, 1> result;
	result.setZero();

	// create quadrature for line and surface integrals
	typedef NiHu::gaussian_quadrature<NiHu::line_domain> line_quad_t;
	typedef NiHu::gaussian_quadrature<domain_t> surf_quad_t;
	unsigned order = 40;
	line_quad_t line_quad(order);
	surf_quad_t surf_quad(order);

	// contour integrals
	for (unsigned i = 0; i < domain_t::num_corners; ++i)
	{
		xi_t xi1 = domain_t::get_corner(i);
		xi_t xi2 = domain_t::get_corner((i+1)%domain_t::num_corners);
		xi_t dxi_eta = (xi2 - xi1) / 2.0;

		for (auto q : line_quad)
		{
			double eta = q.get_xi()(0);
			double w = q.get_w();

			xi_t xi = (1.0 - eta) / 2.0 * xi1 + (1.0 + eta) / 2.0 * xi2;
			typename lset_t::shape_t L = lset_t::template eval_shape<0>(xi);
			typename lset_t::dshape_t dL = lset_t::template eval_shape<1>(xi);
			x_t y = elem.get_coords() * L;
			x_t dy = elem.get_coords() * dL * dxi_eta;

			x_t rvec = y - x;
			double r = rvec.norm();
			double G = 1.0 / (4.0  * pi * r);
			x_t gradG = -1.0 / (4.0  * pi * r * r) * rvec.normalized();

			result += (N0 + (gradN0 * (y - y0))) * w * nx.dot(gradG.cross(dy));
			result -= G * w * gradN0 * dy.cross(nx);
		}

	}

	// surface integrals
	double Om = 0.0;
	for (auto q : surf_quad)
	{
		xi_t xi = q.get_xi();
		double w = q.get_w();
		typename lset_t::shape_t L = lset_t::template eval_shape<0>(xi);
		typename lset_t::dshape_t dL = lset_t::template eval_shape<1>(xi);
		x_t y = elem.get_x(xi);
		x_t jvec = elem.get_normal(xi);
		double jac = jvec.norm();
		x_t ny = jvec.normalized();
		NiHu::laplace_3d_DLPt_kernel dlpt_kernel;
		NiHu::laplace_3d_DLP_kernel dlp_kernel;
		double dGnx = dlpt_kernel(x, y, nx);
		double dGny = dlp_kernel(x, y, ny);
		result += dGnx * w * jac * (gradN0 * ny);
		Om -= 4.0 * pi * w * jac * dGny;
	}

	result -= (gradN0 * nx) * Om / (4.0 * pi);

	return result;
}

int main()
{
	using namespace boost::math::double_constants;
	
	typedef NiHu::laplace_3d_HSP_kernel hsp_kernel_t;
	typedef NiHu::quad_1_elem elem_t;
	typedef elem_t::x_t x_t;
	typedef elem_t::xi_t xi_t;
	typedef elem_t::lset_t lset_t;
	typedef NiHu::quad_1_shape_set nset_t;

	elem_t::coords_t coords;
	coords <<
		-1.0, 1.0, 1.0, -1.0,
		-1.0, -1.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	elem_t elem(coords);
	x_t x;
	x << .5, .5, .2;

	x_t nx;
	nx << 0.0, 0.0, 1.0;

	// inverse mapping
	NiHu::inverse_mapping<elem_t> invmap(elem);
	invmap.eval(x, 1e-6, 100);
	xi_t xi0 = invmap.get_result().topRows(2);

	// local Taylor series expansion
	nset_t::shape_t N0 = nset_t::eval_shape<0>(xi0);
	nset_t::dshape_t dN0 = nset_t::eval_shape<1>(xi0);
	Eigen::Matrix<double, 4, 3> dN0_x;
	dN0_x.leftCols(2) = dN0;
	dN0_x.col(2).setZero();

	nset_t::shape_t L0 = lset_t::eval_shape<0>(xi0);
	nset_t::dshape_t dL0 = lset_t::eval_shape<1>(xi0);
	Eigen::Matrix<double, 3, 3> dy0;
	dy0.col(0) = (coords * dL0.col(0)).transpose();
	dy0.col(1) = (coords * dL0.col(1)).transpose();
	dy0.col(2) = dy0.col(0).cross(dy0.col(1));

	Eigen::Matrix<double, 4, 3> gradN0 = dN0_x * dy0.inverse();

	x_t y0 = coords * L0;

	Eigen::Matrix<double, 4, 1> Istokes = stokes_hsp_integral<elem_t, nset_t>(
		elem, x, nx, y0, N0, gradN0);

	std::cout << std::setprecision(10) << Istokes(0, 0) / N0(0, 0) << std::endl;

	std::cout << std::endl;

	std::cout << std::setprecision(10) << Istokes.sum() * (4.0 * pi) << std::endl;

	return 0;
}
