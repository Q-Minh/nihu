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

#include "library/helmholtz_kernel.hpp"
#include <iostream>

template <class kernel>
void tester(kernel_base<kernel> const &k)
{
	line_1_elem::coords_t coords1;
	coords1 <<
		0, 1,
		0, 0;
	line_1_elem::coords_t coords2;
	coords2 <<
		1, 2,
		1, 1;

	line_1_elem elem1(coords1);
	line_1_elem elem2(coords2);
	auto xi1 = line_1_elem::domain_t::get_center();
	auto xi2 = line_1_elem::domain_t::get_center();

	typename kernel_base<kernel>::test_input_t x(elem1, xi1);
	typename kernel_base<kernel>::trial_input_t y(elem2, xi2);
	std::cout << k(x,y) << std::endl;
}

int main(void)
{
	double wave_number(1.0);
	std::cout << "helmholtz G kernel";
	tester(helmholtz_2d_SLP_kernel<double>(wave_number));
	std::cout << "helmholtz H kernel";
	tester(helmholtz_2d_DLP_kernel<double>(wave_number));

	return 0;
}

