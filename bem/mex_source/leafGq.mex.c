/*leafGq
 *
 *
 * Peter Fiala
 * 2009
 */

/* $Make: mex -O -output leafGq leafGq.mex.c fmbem.c vector.c$ */

#include "mex.h"
#include "fmbem.h"

/* ------------------------------------------------------------------------ */
/* The entry point of the mex function. */
/* This function reads Matlab parameters, allocates space for the output */
/* and calls C routines that perform the computations */
void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
    double *r, *qr, *qi, *R, *s, k;
    int32_t *father;
    int nnod, nclus, ns;
    double *Fr, *Fi;

    /* transfer input parameters */
    r = mxGetPr(prhs[0]);		/* source nodes */
    nnod = mxGetN(prhs[0]);
    qr = mxGetPr(prhs[1]);		/* excitation pressure derivative */
    qi = mxGetPi(prhs[1]);
    father = (int32_t *)mxGetData(prhs[2]);	/* father cluster indices */
    nclus = (int)mxGetScalar(prhs[3]);		/* cluster centres */
    s = mxGetPr(prhs[4]);		/* directions (unit sphere) */
    ns = mxGetN(prhs[4]);
    k = mxGetScalar(prhs[5]);	/* wave number */

    /* Allocate output parameters */
    plhs[0] = mxCreateDoubleMatrix(ns, nclus, mxCOMPLEX);
    Fr = mxGetPr(plhs[0]);		/* Far field signature */
    Fi = mxGetPi(plhs[0]);

    /* call C subroutine */
    leafGq(nnod, r, qr, qi, father, nclus, ns, s, k, Fr, Fi);

    return;
}
