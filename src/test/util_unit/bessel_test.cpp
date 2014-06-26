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
#include "util/math_functions.hpp"

int main(void)
{
	std::cout << bessel::K<0>(std::complex<double>(0.)) << std::endl;
	std::cout << bessel::K<0>(std::complex<double>(1.)) << std::endl;
	std::cout << bessel::K<0>(std::complex<double>(2.)) << std::endl;

	std::cout << bessel::K<0>(std::complex<double>(0., 1.)) << std::endl;
	std::cout << bessel::K<0>(std::complex<double>(1., 1.)) << std::endl;
	std::cout << bessel::K<0>(std::complex<double>(2., 1.)) << std::endl;

	std::cout << bessel::K<0>(std::complex<double>(0., -1.)) << std::endl;
	std::cout << bessel::K<0>(std::complex<double>(1., -1.)) << std::endl;
	std::cout << bessel::K<0>(std::complex<double>(2., -1.)) << std::endl;
}

