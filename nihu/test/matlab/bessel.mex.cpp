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

#include "nihu/util/mex_matrix.hpp"
#include "nihu/util/math_functions.hpp"
#include <iostream>
#include <stdexcept>

typedef NiHu::mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	if (nlhs != 10 || nrhs != 1)
		mexErrMsgIdAndTxt("bessel:arguments",
			"bad number of input (%d) or output (%d) arguments", nlhs, nrhs);

	cMatrix z(rhs[0]);
	cMatrix J0(z.rows(), z.cols(), lhs[0]), J1(z.rows(), z.cols(), lhs[1]);
	cMatrix Y0(z.rows(), z.cols(), lhs[2]), Y1(z.rows(), z.cols(), lhs[3]);
	cMatrix H01(z.rows(), z.cols(), lhs[4]), H11(z.rows(), z.cols(), lhs[5]);
	cMatrix H02(z.rows(), z.cols(), lhs[6]), H12(z.rows(), z.cols(), lhs[7]);
	cMatrix K0(z.rows(), z.cols(), lhs[8]), K1(z.rows(), z.cols(), lhs[9]);

	for (unsigned i = 0; i < z.rows(); ++i)
		for (unsigned j = 0; j < z.cols(); ++j)
		{
			J0(i,j) = NiHu::bessel::J<0, std::complex<double> >(z(i,j));
			J1(i,j) = NiHu::bessel::J<1, std::complex<double> >(z(i,j));
			Y0(i,j) = NiHu::bessel::Y<0>(z(i,j));
			Y1(i,j) = NiHu::bessel::Y<1>(z(i,j));
			H01(i,j) = NiHu::bessel::H<0, 1, std::complex<double> >(z(i,j));
			H11(i,j) = NiHu::bessel::H<1, 1, std::complex<double> >(z(i,j));
			H02(i,j) = NiHu::bessel::H<0, 2, std::complex<double> >(z(i,j));
			H12(i,j) = NiHu::bessel::H<1, 2, std::complex<double> >(z(i,j));
			K0(i,j) = NiHu::bessel::K<0, std::complex<double> >(z(i,j));
			K1(i,j) = NiHu::bessel::K<1, std::complex<double> >(z(i,j));
		}
}

