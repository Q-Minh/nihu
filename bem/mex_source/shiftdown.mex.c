/*SHIFTDOWN
 *
 *
 * Peter Fiala
 * 2009
 */

/* $Make: mex -O -output shiftdown shiftdown.mex.c fmbem.c vector.c$ */

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
    double *r, *qr, *qi, *father, *R, *s, k;
    int nnod, nclus, ns;
    double *Nr, *Ni;

    /* transfer input parameters */
    r = mxGetPr(prhs[0]);		/* source nodes */
    nnod = mxGetN(prhs[0]);
    qr = mxGetPr(prhs[1]);		/* excitation */
    qi = mxGetPi(prhs[1]);
    nclus = (int)mxGetScalar(prhs[1]);
    father = mxGetPr(prhs[2]);	/* father cluster indices */
    s = mxGetPr(prhs[3]);		/* directions (unit sphere) */
    ns = mxGetN(prhs[3]);
    k = mxGetScalar(prhs[4]);	/* wave number */

    /* output vector */
    plhs[0] = mxCreateDoubleMatrix(ns, nnod, mxCOMPLEX);
    Nr = mxGetPr(plhs[0]);		/* Near field signature */
    Ni = mxGetPi(plhs[0]);

    /* call C subroutine */
    downward(nnod, r, qr, qi, father, ns, s, k, Nr, Ni);

    return;
}
