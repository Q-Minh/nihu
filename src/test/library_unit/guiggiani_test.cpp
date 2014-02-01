#include <iostream>
#include <iomanip>

#include "library/guiggiani_1992.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"
#include "library/quad_28_elem.hpp"

template <class elem_t>
double laplace_planar_analytic(typename elem_t::coords_t const &coords, typename elem_t::x_t const &x0)
{
	typedef typename elem_t::x_t x_t;
	elem_t elem(coords);
	double res = 0.0;
	unsigned const N = elem_t::domain_t::num_corners;
	for (unsigned i = 0; i < N; ++i)
	{
		x_t c1 = elem.get_coords().col(i);
		x_t c2 = elem.get_coords().col((i + 1) % N);
		x_t d1 = c1 - x0, d2 = c2 - x0;
		x_t side = (c2 - c1).normalized();
		x_t d = d1 - side * side.dot(d1);

		double phi1 = std::asin(d1.normalized().cross(d.normalized())(2));
		double phi2 = std::asin(d.normalized().cross(d2.normalized())(2));

		res += (std::sin(phi2) + std::sin(phi1)) / d.norm();
	}
	return -res / (4.0 * M_PI);
}


template <class Kernel, class Elem, unsigned order>
struct tester
{
	static typename Kernel::scalar_t
	eval(Kernel const &kernel, typename Elem::coords_t const &coords, typename Elem::xi_t const &xi0)
	{
		Elem elem(coords);
		guiggiani<
			field_view<Elem, field_option::constant>,
			Kernel,
			2 * order - 1
		> gui(elem, kernel);
		Eigen::Matrix<typename Kernel::scalar_t, 1, 1> result;
		result.setZero();
		gui.integrate(result, xi0, elem.get_normal(xi0).normalized());
		return result(0, 0);
	}
};


#define TEST \
for (unsigned i = 0; i < sizeof(xi0) / sizeof(xi0[0]); ++i) \
{\
	quad_1_elem elem(coords); \
	I0[i] = laplace_planar_analytic<quad_1_elem>(coords, elem.get_x(xi0[i])); \
	I[i][0] = tester<laplace_3d_HSP_kernel, quad_1_elem, 3>::eval(kernel, coords, xi0[i]);\
	I[i][1] = tester<laplace_3d_HSP_kernel, quad_1_elem, 5>::eval(kernel, coords, xi0[i]);\
	I[i][2] = tester<laplace_3d_HSP_kernel, quad_1_elem, 7>::eval(kernel, coords, xi0[i]);\
	I[i][3] = tester<laplace_3d_HSP_kernel, quad_1_elem, 9>::eval(kernel, coords, xi0[i]);\
}\
std::cout << "xi0:";\
for (unsigned i = 0; i < sizeof(xi0) / sizeof(xi0[0]); ++i)\
	std::cout << "\t(" << xi0[i](0) << "," << xi0[i](1) << ")";\
std::cout << std::endl;\
for (unsigned j = 0; j < 4; ++j)\
{\
	std::cout << "Err:";\
	for (unsigned i = 0; i < sizeof(xi0) / sizeof(xi0[0]); ++i)\
		std::cout << "\t" << std::log10(std::abs(I[i][j] / I0[i] - 1.0));\
	std::cout << std::endl;\
}\



void test_laplace_3d_plane_linear(void)
{
	laplace_3d_HSP_kernel kernel;
	quad_1_elem::coords_t coords;
	quad_1_elem::xi_t xi0[3];
	xi0[0] << 0.0, 0.0;
	xi0[1] << 0.57, 0.0;
	xi0[2] << 0.9, 0.923;

	double I0[3], I[3][4];

	std::cout << std::setprecision(4);

	std::cout << "\nLinear quad element:\n===================\n";

	std::cout << "\nSquare:\n----------------\n";

	coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	TEST

	std::cout << "\nRectangle:\n----------------\n";

	coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	coords.row(0) *= 2;

	TEST

	std::cout << "\nParallelogram:\n----------------\n";

	coords <<
		0.0, 1.0, 3.0, 2.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	TEST

	std::cout << "\nSlightly distorted:\n----------------\n";

	coords <<
		0.0, 0.4, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	TEST

	std::cout << "\nHighly distorted:\n----------------\n";

	coords <<
		0.0, 0.1, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	TEST

}

#undef TEST

#define TEST \
for (unsigned i = 0; i < sizeof(xi0) / sizeof(xi0[0]); ++i) \
{\
	quad_2_elem elem(coords); \
	I0[i] = laplace_planar_analytic<quad_1_elem>(coords0, elem.get_x(xi0[i])); \
	I[i][0] = tester<laplace_3d_HSP_kernel, quad_2_elem, 3>::eval(kernel, coords, xi0[i]); \
	I[i][1] = tester<laplace_3d_HSP_kernel, quad_2_elem, 5>::eval(kernel, coords, xi0[i]); \
	I[i][2] = tester<laplace_3d_HSP_kernel, quad_2_elem, 7>::eval(kernel, coords, xi0[i]); \
	I[i][3] = tester<laplace_3d_HSP_kernel, quad_2_elem, 9>::eval(kernel, coords, xi0[i]); \
	}\
	std::cout << "xi0:"; \
for (unsigned i = 0; i < sizeof(xi0) / sizeof(xi0[0]); ++i)\
	std::cout << "\t(" << xi0[i](0) << "," << xi0[i](1) << ")"; \
	std::cout << std::endl; \
for (unsigned j = 0; j < 4; ++j)\
{\
	std::cout << "Err:"; \
for (unsigned i = 0; i < sizeof(xi0) / sizeof(xi0[0]); ++i)\
	std::cout << "\t" << std::log10(std::abs(I[i][j] / I0[i] - 1.0)); \
	std::cout << std::endl; \
	}\



void test_laplace_3d_plane_quadratic(void)
{
	laplace_3d_HSP_kernel kernel;
	quad_2_elem::coords_t coords;
	quad_1_elem::coords_t coords0;
	quad_2_elem::xi_t xi0[3];
	xi0[0] << 0.0, 0.0;
	xi0[1] << 0.57, 0.0;
	xi0[2] << 0.9, 0.923;

	double I0[3], I[3][4];

	std::cout << std::setprecision(4);

	std::cout << "\nQuadratic quad element:\n===================\n";

	std::cout << "\nSquare:\n----------------\n";

	coords <<
		0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0, 0.5,
		0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.5,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	coords0 <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	TEST

	std::cout << "\nRectangle:\n----------------\n";

	coords.row(0) *= 2;
	coords0.row(0) *= 2;

	TEST

	std::cout << "\nParallelogram:\n----------------\n";

	coords <<
		0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0, 0.5,
		0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.5,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	coords0 <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
	coords.row(0) = 2 * coords.row(0) + coords.row(1);
	coords0.row(0) = 2 * coords0.row(0) + coords0.row(1);

	TEST

	std::cout << "\nDistorted square:\n----------------\n";

	coords <<
		0.0, 0.4, 1.0, 1.0, 1.0, 0.6, 0.0, 0.0, 0.5,
		0.0, 0.0, 0.0, 0.7, 1.0, 1.0, 1.0, 0.3, 0.35,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	coords0 <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	TEST

}


void test_rong(void)
{
	typedef double wave_number_t;
	typedef helmholtz_3d_HSP_kernel<wave_number_t> kernel_t;
	kernel_t kernel(0.0);

	typedef tria_2_elem elem_t;
	typedef field_view<elem_t, field_option::constant> field_t;

	typedef elem_t::coords_t coords_t;

	coords_t coords;
	coords <<
		1.0, 1.0, 1.0, std::cos(M_PI / 8.0), std::cos(M_PI / 4.0), std::cos(M_PI / 8.0),
		-0.5, 0.0, 0.5, 0.25, 0.0, -0.25,
		0.0, 0.0, 0.0, std::sin(M_PI / 8.0), std::sin(M_PI / 4.0), std::sin(M_PI / 8.0);

	elem_t elem(coords);

	elem_t::xi_t xi0(.3, .3);

	Eigen::Matrix<std::complex<double>, 1, 1> I5, I8, I10, I12, I20;

	guiggiani<field_t, kernel_t, 2 * 20 - 1> gui20(elem, kernel);
	I20.setZero();
	gui20.integrate(I20, xi0, elem.get_normal(xi0).normalized());

	guiggiani<field_t, kernel_t, 2 * 5 - 1> gui5(elem, kernel);
	I5.setZero();
	gui5.integrate(I5, xi0, elem.get_normal(xi0).normalized());
	std::cout << I5 << ' ' << std::log10(std::abs(I5(0, 0) / I20(0, 0) - 1.0)) << std::endl;

	guiggiani<field_t, kernel_t, 2 * 8 - 1> gui8(elem, kernel);
	I8.setZero();
	gui8.integrate(I8, xi0, elem.get_normal(xi0).normalized());
	std::cout << I8 << ' ' << std::log10(std::abs(I8(0, 0) / I20(0, 0) - 1.0)) << std::endl;

	guiggiani<field_t, kernel_t, 2 * 10 - 1> gui10(elem, kernel);
	I10.setZero();
	gui10.integrate(I10, xi0, elem.get_normal(xi0).normalized());
	std::cout << I10 << ' ' << std::log10(std::abs(I10(0, 0) / I20(0, 0) - 1.0)) << std::endl;

	guiggiani<field_t, kernel_t, 2 * 12 - 1> gui12(elem, kernel);
	I12.setZero();
	gui12.integrate(I12, xi0, elem.get_normal(xi0).normalized());
	std::cout << I12 << ' ' << std::log10(std::abs(I12(0, 0) / I20(0, 0) - 1.0)) << std::endl;

}


void test_guiggiani_92_curved(void)
{
	typedef double wave_number_t;
	typedef laplace_3d_HSP_kernel kernel_t;
	kernel_t kernel;

	typedef quad_2_elem elem_t;
	typedef field_view<elem_t, field_option::constant> field_t;

	typedef elem_t::coords_t coords_t;

	coords_t coords;
	double b = std::sqrt(2.0) / 2.0;
	coords <<
		1, 1, 1, b, 0, 0, 0, b, b,
		0, 1, 2, 2, 2, 1, 0, 0, 1,
		0, 0, 0, b, 1, 1, 1, b, b;

	elem_t elem(coords);
	guiggiani<field_t, kernel_t, 7> gui(elem, kernel);

	elem_t::xi_t xi0;
	Eigen::Matrix<double, 1, 1> I;

	I.setZero();
	xi0 << 0.0, 0.0;
	gui.integrate(I, xi0, elem.get_normal(xi0).normalized());
	std::cout << I << std::endl;

	//I.setZero();
	//xi0 << 0.66, 0.0;
	//gui.integral(I, xi0);
	//std::cout << I << std::endl;

	//I.setZero();
	//xi0 << 0.66, 0.66;
	//gui.integral(I, xi0);
	//std::cout << I << std::endl;
}


void test_helmholtz_3d(void)
{
	std::cout << "distortion test of helmholtz 3D kernel with quadratic triangle:\n";

	typedef tria_1_elem elem0_t;
	typedef tria_2_elem elem_t;
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef trial_field_t test_field_t;
	typedef helmholtz_3d_HSP_kernel<double> kernel_t;

	elem0_t::coords_t coords0;
	coords0 <<
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 0.0, 0.0;

	elem_t::coords_t coords;
	coords <<
		0.0, 0.5, 1.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.5, 1.0, 0.5,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	coords0.row(1) *= 2.0;
	coords.row(1) *= 2.0;

	elem0_t elem0(coords0);
	elem_t elem(coords);

	double wave_num = .1;

	kernel_t kernel(wave_num);

	guiggiani<trial_field_t, kernel_t, 5> gui(elem, kernel);

	Eigen::Matrix<std::complex<double>, 1, 1> I;
	I.setZero();
	auto xi0 = elem_t::domain_t::get_center();
	gui.integrate(I, xi0, elem.get_normal(xi0).normalized());
	std::complex<double> I0 = helmholtz_3d_HSP_collocation_constant_triangle::eval(
		elem0,
		elem.get_center(),
		wave_num);

	std::cout << "I:\t" << I << std::endl;
	std::cout << "Ianal:\t" << I0 << std::endl;
	std::cout << "log10 error:\t" << std::log10(std::abs((I / I0).norm() - 1.0)) << std::endl;
}

int main(void)
{
	test_laplace_3d_plane_linear();
	test_laplace_3d_plane_quadratic();

	//	test_laplace_3d_quadratic();
	//	test_helmholtz_3d();

	test_guiggiani_92_curved();
	test_rong();

	return 0;
}
