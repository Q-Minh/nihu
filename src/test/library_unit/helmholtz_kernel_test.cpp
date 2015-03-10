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

#include "library/helmholtz_kernel.hpp"
#include "library/lib_element.hpp"
#include <iostream>

void comparison()
{
	std::complex<double> const I(0., 1.);
	
	double x(1.0), y(1.0), r = std::sqrt(x*x+y*y);
	std::complex<double> s(1., -1.);
	std::complex<double>k(-I*s);
	
	std::complex<double> krH0((k*r)*bessel::H<0,2>(k*r));
	std::complex<double> H1(bessel::H<1,2>(k*r));
	
	std::complex<double> rsK0(r*s*bessel::K<0>(r*s));
	std::complex<double> K1(bessel::K<1>(r*s));

	std::complex<double> Ker1 = I*k/r * (x/r*x/r*(krH0-2.*H1) + H1) / 4.;
	std::complex<double> Ker2 = s/r * (x/r*x/r*(2.*K1+rsK0) - K1) / (2.*M_PI);
	
	std::cout << Ker1 << '\t' << Ker2 << std::endl;
}

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

	typename kernel_traits<kernel>::test_input_t x(elem1, xi1);
	typename kernel_traits<kernel>::trial_input_t y(elem2, xi2);

	std::cout << k(x,y) << std::endl;
}

int main(void)
{
	std::complex<double> wave_number(-1., -1.);

	std::cout << "helmholtz G kernel";
	tester(helmholtz_2d_SLP_kernel<std::complex<double> >(wave_number));

	std::cout << "helmholtz H kernel";
	tester(helmholtz_2d_DLP_kernel<std::complex<double> >(wave_number));
	
	std::cout << "helmholtz double kernel";
	tester(helmholtz_2d_double_kernel<std::complex<double>, 0, 0>(wave_number));
	
	comparison();

	return 0;
}

