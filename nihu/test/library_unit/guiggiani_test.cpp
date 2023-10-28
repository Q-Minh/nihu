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

#include "nihu/library/guiggiani_1992.hpp"
#include "nihu/library/laplace_nearly_singular_integrals.hpp"
#include "nihu/library/laplace_singular_integrals.hpp"
#include "nihu/library/helmholtz_nearly_singular_integrals.hpp"
#include "nihu/library/helmholtz_singular_integrals.hpp"
#include "nihu/library/lib_element.hpp"

template <class kernel_t>
struct planar_analytic;

template <class wavenumber_t>
struct planar_analytic<NiHu::helmholtz_3d_HSP_kernel<wavenumber_t> >
{
	template <class elem_t>
	static typename NiHu::helmholtz_3d_HSP_kernel<wavenumber_t>::result_t
	eval(NiHu::helmholtz_3d_HSP_kernel<wavenumber_t> const &kernel, elem_t const &elem, typename elem_t::x_t const &x0)
	{
		return NiHu::helmholtz_3d_HSP_collocation_constant_plane<40>::eval(elem, x0, kernel.get_wave_number());
	}
};

template <>
struct planar_analytic<NiHu::laplace_3d_HSP_kernel>
{
	template <class elem_t>
	static typename NiHu::laplace_3d_HSP_kernel::result_t
		eval(NiHu::kernel_base<NiHu::laplace_3d_HSP_kernel> const &kernel, elem_t const &elem, typename elem_t::x_t const &x0)
	{
		return NiHu::laplace_3d_HSP_collocation_constant_plane::eval(elem, x0);
	}
};


template <class Kernel, class Elem, unsigned order>
struct tester
{
	static typename Kernel::result_t
	eval(NiHu::kernel_base<Kernel> const &kernel, typename Elem::coords_t const &coords, typename Elem::xi_t const &xi0)
	{
		Elem elem(coords);
		NiHu::guiggiani<
			NiHu::field_view<Elem, NiHu::field_option::constant>,
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
	NiHu::quad_1_elem elem(coords); \
	NiHu::quad_1_elem elem0 = elem; \
	I0[i] = planar_analytic<kernel_t>::eval(kernel, elem0, elem0.get_x(xi0[i])); \
	I[i][0] = tester<kernel_t, NiHu::quad_1_elem, 3>::eval(kernel, coords, xi0[i]);\
	I[i][1] = tester<kernel_t, NiHu::quad_1_elem, 5>::eval(kernel, coords, xi0[i]);\
	I[i][2] = tester<kernel_t, NiHu::quad_1_elem, 7>::eval(kernel, coords, xi0[i]);\
	I[i][3] = tester<kernel_t, NiHu::quad_1_elem, 9>::eval(kernel, coords, xi0[i]);\
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
	// corner coordinates of the element
	NiHu::quad_1_elem::coords_t coords;
	
	// intrinsic coordinates of the singular points
	NiHu::quad_1_elem::xi_t xi0[3];
	xi0[0] << 0.0, 0.0;
	xi0[1] << 0.57, 0.0;
	xi0[2] << 0.9, 0.923;

	// I0 is the analytic solution, I contains the numerical solutions
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
	NiHu::quad_2_elem elem(coords); \
	NiHu::quad_1_elem elem0(coords0); \
	I0[i] = planar_analytic<kernel_t>::eval(kernel, elem0, elem.get_x(xi0[i])); \
	I[i][0] = tester<kernel_t, NiHu::quad_2_elem, 3>::eval(kernel, coords, xi0[i]); \
	I[i][1] = tester<kernel_t, NiHu::quad_2_elem, 5>::eval(kernel, coords, xi0[i]); \
	I[i][2] = tester<kernel_t, NiHu::quad_2_elem, 7>::eval(kernel, coords, xi0[i]); \
	I[i][3] = tester<kernel_t, NiHu::quad_2_elem, 9>::eval(kernel, coords, xi0[i]); \
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
	NiHu::quad_2_elem::coords_t coords;
	NiHu::quad_1_elem::coords_t coords0;
	NiHu::quad_2_elem::xi_t xi0[3];
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


int main(void)
{
	std::cout << "\n\nTesting with Laplace kernel\n";

	NiHu::laplace_3d_HSP_kernel l_kernel;
	test_plane_linear(l_kernel);
	test_plane_quadratic(l_kernel);

	std::cout << "\n\nTesting with Helmholtz kernel\n";

	NiHu::helmholtz_3d_HSP_kernel<double> h_kernel(.2);
	test_plane_linear(h_kernel);
	test_plane_quadratic(h_kernel);

	return 0;
}

