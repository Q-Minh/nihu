/*RECOVER
 *
 *
 * Peter Fiala
 * 2009
 */

/* $Make: mex -O -output recover recover.mex.c fmbem.c vector.c$ */

#include "fmbem.h"
#include "mex.h"

/* ------------------------------------------------------------------------ */
/* The entry point of the mex function. */
/* This function reads Matlab parameters, allocates space for the output */
/* and calls C routines that perform the computations */
void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
    double *Nr, *Ni, *father, *r, *s, *w, k;
    int nnod, ns;
    double *pr, *pi;

    /* transfer input parameters */
    r = mxGetPr(prhs[0]);		/* receiver nodes */
    nnod = mxGetN(prhs[0]);
    Nr = mxGetPr(prhs[1]);		/* Near field signature */
    Ni = mxGetPi(prhs[1]);
    father = mxGetPr(prhs[2]);	/* receiver father cluster indices */
    s = mxGetPr(prhs[3]);		/* directions (unit sphere) */
    ns = mxGetN(prhs[3]);
    w = mxGetPr(prhs[4]);		/* Gauss weights over unit sphere */
    k = mxGetScalar(prhs[5]);	/* wave number */

    /* output vector */
    plhs[0] = mxCreateDoubleMatrix(nnod, 1, mxCOMPLEX);
    pr = mxGetPr(plhs[0]);		/* response pressures */
    pi = mxGetPi(plhs[0]);

    /* call C subroutine */
    recover(nnod, r, Nr, Ni, father, ns, s, w, k, pr, pi);

    return;
}
