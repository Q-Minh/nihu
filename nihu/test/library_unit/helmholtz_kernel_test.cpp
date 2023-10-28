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

#include <boost/math/constants/constants.hpp>

#include "nihu/library/helmholtz_kernel.hpp"
#include "nihu/library/lib_element.hpp"
#include <iostream>


void comparison()
{
	using namespace boost::math::double_constants;

	std::complex<double> const I(0., 1.);
	
	double x(1.0), y(1.0), r = std::sqrt(x*x+y*y);
	std::complex<double> s(1., -1.);
	std::complex<double>k(-I*s);
	
	std::complex<double> krH0((k*r)*NiHu::bessel::H<0,2>(k*r));
	std::complex<double> H1(NiHu::bessel::H<1,2>(k*r));
	
	std::complex<double> rsK0(r*s*NiHu::bessel::K<0>(r*s));
	std::complex<double> K1(NiHu::bessel::K<1>(r*s));

	std::complex<double> Ker1 = I*k/r * (x/r*x/r*(krH0-2.*H1) + H1) / 4.;
	std::complex<double> Ker2 = s/r * (x/r*x/r*(2.*K1+rsK0) - K1) / two_pi;
	
	std::cout << Ker1 << '\t' << Ker2 << std::endl;
}


template <class kernel>
void tester(NiHu::kernel_base<kernel> const &k)
{
	NiHu::line_1_elem::coords_t coords1;
	coords1 <<
		0, 1,
		0, 0;

	NiHu::line_1_elem::coords_t coords2;
	coords2 <<
		1, 2,
		1, 1;

	NiHu::line_1_elem elem1(coords1);
	NiHu::line_1_elem elem2(coords2);
	auto xi1 = NiHu::line_1_elem::domain_t::get_center();
	auto xi2 = NiHu::line_1_elem::domain_t::get_center();

	typename NiHu::kernel_traits<kernel>::test_input_t x(elem1, xi1);
	typename NiHu::kernel_traits<kernel>::trial_input_t y(elem2, xi2);

	std::cout << k(x,y) << std::endl;
}

int main(void)
{
	std::complex<double> wave_number(-1., -1.);

	std::cout << "helmholtz G kernel";
	tester(NiHu::helmholtz_2d_SLP_kernel<std::complex<double> >(wave_number));

	std::cout << "helmholtz H kernel";
	tester(NiHu::helmholtz_2d_DLP_kernel<std::complex<double> >(wave_number));
	
	comparison();

	return 0;
}

