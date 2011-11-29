/* $Make: mex -O -output mex_test mex_test.mex.c$ */

#include "mex.h"

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
    /* input parameters */
    double *in, *out;

    in = mxGetPr(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    out = mxGetPr(plhs[0]);
    *out = 2.0 * *in;
    return;
}
