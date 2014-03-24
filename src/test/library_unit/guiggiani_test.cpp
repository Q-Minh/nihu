// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>

#include "library/guiggiani_1992.hpp"
#include "library/laplace_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"
#include "library/lib_element.hpp"

template <class kernel_t>
struct planar_analytic;

template <class wavenumber_t>
struct planar_analytic<helmholtz_3d_HSP_kernel<wavenumber_t> >
{
	template <class elem_t>
	static typename helmholtz_3d_HSP_kernel<wavenumber_t>::result_t
	eval(helmholtz_3d_HSP_kernel<wavenumber_t> const &kernel, elem_t const &elem, typename elem_t::x_t const &x0)
	{
		return helmholtz_3d_HSP_collocation_constant_plane<40>::eval(elem, x0, kernel.get_data().get_wave_number());
	}
};

template <>
struct planar_analytic<laplace_3d_HSP_kernel>
{
	template <class elem_t>
	static typename laplace_3d_HSP_kernel::result_t
		eval(kernel_base<laplace_3d_HSP_kernel> const &kernel, elem_t const &elem, typename elem_t::x_t const &x0)
	{
		return laplace_3d_HSP_collocation_constant_plane::eval(elem, x0);
	}
};


template <class Kernel, class Elem, unsigned order>
struct tester
{
	static typename Kernel::result_t
	eval(kernel_base<Kernel> const &kernel, typename Elem::coords_t const &coords, typename Elem::xi_t const &xi0)
	{
		Elem elem(coords);
		guiggiani<
			field_view<Elem, field_option::constant>,
			Kernel,
			2 * order - 1
		> gui(elem, kernel);
		Eigen::Matrix<typename Kernel::result_t, 1, 1> result;
		result.setZero();
		gui.integrate(result, xi0, elem.get_normal(xi0).normalized());
		return result(0, 0);
	}
};


#define TEST \
	for (unsigned i = 0; i < sizeof(xi0) / sizeof(xi0[0]); ++i) \
{\
	quad_1_elem elem(coords); \
	quad_1_elem elem0 = elem; \
	I0[i] = planar_analytic<kernel_t>::eval(kernel, elem0, elem0.get_x(xi0[i])); \
	I[i][0] = tester<kernel_t, quad_1_elem, 3>::eval(kernel, coords, xi0[i]);\
	I[i][1] = tester<kernel_t, quad_1_elem, 5>::eval(kernel, coords, xi0[i]);\
	I[i][2] = tester<kernel_t, quad_1_elem, 7>::eval(kernel, coords, xi0[i]);\
	I[i][3] = tester<kernel_t, quad_1_elem, 9>::eval(kernel, coords, xi0[i]);\
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


template <class kernel_t>
void test_plane_linear(kernel_t const &kernel)
{
	quad_1_elem::coords_t coords;
	quad_1_elem::xi_t xi0[3];
	xi0[0] << 0.0, 0.0;
	xi0[1] << 0.57, 0.0;
	xi0[2] << 0.9, 0.923;

	typename kernel_t::result_t I0[3], I[3][4];

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
	quad_1_elem elem0(coords0); \
	I0[i] = planar_analytic<kernel_t>::eval(kernel, elem0, elem.get_x(xi0[i])); \
	I[i][0] = tester<kernel_t, quad_2_elem, 3>::eval(kernel, coords, xi0[i]); \
	I[i][1] = tester<kernel_t, quad_2_elem, 5>::eval(kernel, coords, xi0[i]); \
	I[i][2] = tester<kernel_t, quad_2_elem, 7>::eval(kernel, coords, xi0[i]); \
	I[i][3] = tester<kernel_t, quad_2_elem, 9>::eval(kernel, coords, xi0[i]); \
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



template <class kernel_t>
void test_plane_quadratic(kernel_t const &kernel)
{
	quad_2_elem::coords_t coords;
	quad_1_elem::coords_t coords0;
	quad_2_elem::xi_t xi0[3];
	xi0[0] << 0.0, 0.0;
	xi0[1] << 0.57, 0.0;
	xi0[2] << 0.9, 0.923;

	typename kernel_t::result_t I0[3], I[3][4];

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

//
//void test_rong(void)
//{
//	typedef double wave_number_t;
//	typedef helmholtz_3d_HSP_kernel<wave_number_t> kernel_t;
//	kernel_t kernel(0.0);
//
//	typedef tria_2_elem elem_t;
//	typedef field_view<elem_t, field_option::constant> field_t;
//
//	typedef elem_t::coords_t coords_t;
//
//	coords_t coords;
//	coords <<
//		1.0, 1.0, 1.0, std::cos(M_PI / 8.0), std::cos(M_PI / 4.0), std::cos(M_PI / 8.0),
//		-0.5, 0.0, 0.5, 0.25, 0.0, -0.25,
//		0.0, 0.0, 0.0, std::sin(M_PI / 8.0), std::sin(M_PI / 4.0), std::sin(M_PI / 8.0);
//
//	elem_t elem(coords);
//
//	elem_t::xi_t xi0(.3, .3);
//
//	Eigen::Matrix<std::complex<double>, 1, 1> I5, I8, I10, I12, I20;
//
//	guiggiani<field_t, kernel_t, 2 * 20 - 1> gui20(elem, kernel);
//	I20.setZero();
//	gui20.integrate(I20, xi0, elem.get_normal(xi0).normalized());
//
//	guiggiani<field_t, kernel_t, 2 * 5 - 1> gui5(elem, kernel);
//	I5.setZero();
//	gui5.integrate(I5, xi0, elem.get_normal(xi0).normalized());
//	std::cout << I5 << ' ' << std::log10(std::abs(I5(0, 0) / I20(0, 0) - 1.0)) << std::endl;
//
//	guiggiani<field_t, kernel_t, 2 * 8 - 1> gui8(elem, kernel);
//	I8.setZero();
//	gui8.integrate(I8, xi0, elem.get_normal(xi0).normalized());
//	std::cout << I8 << ' ' << std::log10(std::abs(I8(0, 0) / I20(0, 0) - 1.0)) << std::endl;
//
//	guiggiani<field_t, kernel_t, 2 * 10 - 1> gui10(elem, kernel);
//	I10.setZero();
//	gui10.integrate(I10, xi0, elem.get_normal(xi0).normalized());
//	std::cout << I10 << ' ' << std::log10(std::abs(I10(0, 0) / I20(0, 0) - 1.0)) << std::endl;
//
//	guiggiani<field_t, kernel_t, 2 * 12 - 1> gui12(elem, kernel);
//	I12.setZero();
//	gui12.integrate(I12, xi0, elem.get_normal(xi0).normalized());
//	std::cout << I12 << ' ' << std::log10(std::abs(I12(0, 0) / I20(0, 0) - 1.0)) << std::endl;
//
//}
//
//
//void test_guiggiani_92_curved(void)
//{
//	typedef laplace_3d_HSP_kernel kernel_t;
//	kernel_t kernel;
//
//	typedef quad_2_elem elem_t;
//	typedef field_view<elem_t, field_option::constant> field_t;
//
//	typedef elem_t::coords_t coords_t;
//
//	coords_t coords;
//	double b = std::sqrt(2.0) / 2.0;
//	coords <<
//		1, 1, 1, b, 0, 0, 0, b, b,
//		0, 1, 2, 2, 2, 1, 0, 0, 1,
//		0, 0, 0, b, 1, 1, 1, b, b;
//
//	elem_t elem(coords);
//	guiggiani<field_t, kernel_t, 7> gui(elem, kernel);
//
//	elem_t::xi_t xi0;
//	Eigen::Matrix<double, 1, 1> I;
//
//	I.setZero();
//	xi0 << 0.0, 0.0;
//	gui.integrate(I, xi0, elem.get_normal(xi0).normalized());
//	std::cout << I << std::endl;
//
//	I.setZero();
//	xi0 << 0.66, 0.0;
//	gui.integral(I, xi0);
//	std::cout << I << std::endl;
//
//	I.setZero();
//	xi0 << 0.66, 0.66;
//	gui.integral(I, xi0);
//	std::cout << I << std::endl;
//}


int main(void)
{
	std::cout << "\n\nTesting with Laplace kernel\n";

	laplace_3d_HSP_kernel l_kernel;
	test_plane_linear(l_kernel);
	test_plane_quadratic(l_kernel);

	std::cout << "\n\nTesting with Helmholtz kernel\n";

	helmholtz_3d_HSP_kernel<double> h_kernel(.2);
	test_plane_linear(h_kernel);
	test_plane_quadratic(h_kernel);

//	test_guiggiani_92_curved();
//	test_rong();

	return 0;
}
