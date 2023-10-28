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

#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/lib_element.hpp"
#include <iostream>
#include <chrono>

unsigned const N = 1000000;

template <class kernel>
void tester(NiHu::kernel_base<kernel> const &k)
{
	NiHu::quad_1_elem::coords_t coords;
	coords <<
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;

	NiHu::quad_1_elem elem(coords);
	auto xi1 = NiHu::quad_1_elem::domain_t::get_corner(0);
	auto xi2 = NiHu::quad_1_elem::domain_t::get_corner(1);

	auto start = std::chrono::steady_clock::now();
	for (unsigned n = 0; n < N; ++n)
	{
		typename NiHu::kernel_base<kernel>::test_input_t x(elem, xi1);
		typename NiHu::kernel_base<kernel>::trial_input_t y(elem, xi2);
		k(x,y);
	}
	auto stop = std::chrono::steady_clock::now();
	auto diff = stop - start;
	std::cout << std::chrono::duration<double, std::micro>(diff).count() << " usec" << std::endl ;
}

int main(void)
{
	std::cout << N << " evaluations: " << std::endl;

	std::cout << "laplace G kernel with wall output:      | ";
	tester(NiHu::laplace_3d_SLP_kernel());

	std::cout << "laplace H kernel with wall output:      | ";
	tester(NiHu::laplace_3d_DLP_kernel());

	std::cout << "laplace H kernel with wall output:      | ";
	tester(NiHu::laplace_3d_DLPt_kernel());

	return 0;
}

