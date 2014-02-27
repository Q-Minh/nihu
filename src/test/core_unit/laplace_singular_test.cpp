// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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

#include <core/weighted_residual.hpp>
#include <library/laplace_kernel.hpp>
#include "library/lib_element.hpp"
// #include <library/laplace_singular_integrals.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

typedef laplace_2d_SLP_kernel kernel_t;
typedef line_1_elem elem_t;
typedef elem_t::coords_t coords_t;

auto kernel = kernel_t();

/*
void constant_collocation_tester(void)
{
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef dirac_field<trial_field_t> test_field_t;
	typedef double_integral<kernel_t, test_field_t, trial_field_t> integral_t;

	coords_t n1;
	n1 <<
		0.0, 1.0,
		0.0, 0.0;
	elem_t elem(n1);

	trial_field_t const &trial_field = create_field_view(elem, field_option::constant());
	test_field_t const &test_field = dirac(trial_field);

	std::cout << integral_t::eval(kernel, test_field, trial_field, std::true_type()) << std::endl;
}
*/

/*
void galerkin_face_tester(void)
{
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef trial_field_t test_field_t;

	coords_t n1;
	n1 <<
		0.0, 1.0,
		0.0, 0.0;
	elem_t elem(n1);

	trial_field_t const &trial_field = create_field_view(elem, field_option::constant());
	test_field_t const &test_field = trial_field;

	typedef double_integral<kernel_t, test_field_t, trial_field_t> integral_t;
	std::cout << integral_t::eval(kernel, test_field, trial_field, std::true_type()) << std::endl;
}
*/

int main(void)
{
/*
	constant_collocation_tester();
*/

	return 0;
}
