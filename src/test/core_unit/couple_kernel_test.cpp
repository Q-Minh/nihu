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

#include <iostream>
#include "core/kernel.hpp"
#include "library/laplace_kernel.hpp"
#include "library/lib_element.hpp"

int main(void)
{
	// create a set of kernels
	auto K = create_couple_kernel(
		laplace_3d_SLP_kernel(),
		laplace_3d_DLP_kernel(),
		laplace_3d_SLP_kernel()
	);

	// create an element
	quad_1_elem::coords_t coords;
	coords <<
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;
	quad_1_elem elem(coords);

	// create test and trial kernel inputs
	auto xi1 = quad_domain::get_center();
	auto xi2 = quad_domain::get_corner(0);
	kernel_traits<decltype(K)>::test_input_t tst_i(elem, xi1);
	kernel_traits<decltype(K)>::trial_input_t trl_i(elem, xi2);

	// evaluate kernel
	auto res = K(tst_i, trl_i);

	// weight kernel by a matrix
	double a = 3.0;
	auto wres = a * res * 2;

	std::cout << wres.get<0>() << std::endl;
	std::cout << wres.get<1>() << std::endl;
	std::cout << wres.get<2>() << std::endl;

	return 0;
}

