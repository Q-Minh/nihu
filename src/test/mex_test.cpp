#include "mex_matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    mex::real_matrix a(prhs[0]);
    mex::complex_matrix b(a.rows(), a.cols(), plhs[0]);
    
    for (size_t i = 0; i < a.rows(); ++i)
        for (size_t j = 0; j < a.cols(); ++j)
            b(i,j) += a(i,j);
}
