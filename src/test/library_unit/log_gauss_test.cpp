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

#include <library/log_gauss_quadrature.hpp>
#include <iostream>

double const D = .2;

double func(double x, unsigned n)
{
	double res = -std::log(x*D);
	for (unsigned k = 0; k < n; ++k)
		res *= x;
	return res;
}

int main(void)
{
	unsigned N = 7;
	auto quadrature = NiHu::log_gauss_impl<double>(N);

	for (unsigned n = 0; n < 2 * N + 5; ++n)
	{
		double res = 0.0;
		for (unsigned i = 0; i < quadrature.rows(); ++i)
			res += func(quadrature(i, 0), n) * quadrature(i, 1);
		double anal = 1.0 / (n + 1) / (n + 1) - std::log(D)/(n+1);
		std::cout << n << '\t' << std::log10(std::abs(res / anal - 1.0)) << std::endl;
	}

	return 0;
}
