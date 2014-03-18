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

#include "util/mex_matrix.hpp"
#include "util/math_functions.hpp"
#include <iostream>
#include <stdexcept>

typedef mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	if (nlhs != 2 || nrhs != 1)
		throw std::invalid_argument("bad number of input or output arguments");

	cMatrix z(rhs[0]);
	cMatrix H0(z.rows(), z.cols(), lhs[0]), H1(z.rows(), z.cols(), lhs[1]);

	for (unsigned i = 0; i < z.rows(); ++i)
		for (unsigned j = 0; j < z.cols(); ++j)
		{
			H0(i,j) = bessel::H<0>(z(i,j));
			H1(i,j) = bessel::H<1>(z(i,j));
		}
}
