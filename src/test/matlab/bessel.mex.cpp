#include "util/mex_matrix.hpp"
#include "util/math_functions.hpp"

typedef mex::complex_matrix<double> cMatrix;

void mexFunction(
	int nlhs, mxArray *lhs[],
	int nrhs, mxArray const *rhs[])
{
	if (nlhs != 1 || nrhs != 1)
		throw ("bad number of input or output arguments");

	cMatrix z(rhs[0]);
	cMatrix H0(z.rows(), z.cols(), lhs[0]);
	for (unsigned i = 0; i < z.rows(); ++i)
		for (unsigned j = 0; j < z.cols(); ++j)
			H0(i,j) = bessel::H<0>(std::complex<double>(z(i,j)));
}
