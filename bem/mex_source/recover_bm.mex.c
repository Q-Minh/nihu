/*RECOVER
 *
 *
 * Peter Fiala
 * 2009
 */

/* $Make: mex -O -output recover_bm recover_bm.mex.c fmbem.c vector.c$ */

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
    double *Nr, *Ni, *r, *n, *s, *w, alphar, alphai, k;
    int32_t *father;
    int nnod, ns;
    double *pr, *pi;

    /* transfer input parameters */
    r = mxGetPr(prhs[0]);		/* receiver nodes */
    nnod = mxGetN(prhs[0]);
	n = mxGetPr(prhs[1]);
    Nr = mxGetPr(prhs[2]);		/* Near field signature */
    Ni = mxGetPi(prhs[2]);
    father = (int32_t *)mxGetData(prhs[3]);	/* father cluster indices */
    s = mxGetPr(prhs[4]);		/* directions (unit sphere) */
    ns = mxGetN(prhs[4]);
    w = mxGetPr(prhs[5]);		/* Gauss weights over unit sphere */
    k = mxGetScalar(prhs[6]);	/* wave number */
	alphar = *(mxGetPr(prhs[7]));/* coupling coefficient */
	alphai = *(mxGetPi(prhs[7]));

    /* output vector */
    plhs[0] = mxCreateDoubleMatrix(nnod, 1, mxCOMPLEX);
    pr = mxGetPr(plhs[0]);		/* response pressures */
    pi = mxGetPi(plhs[0]);

    /* call C subroutine */
    recover_bm(nnod, r, n, Nr, Ni, father, ns, s, w, k, alphar, alphai, pr, pi);

    return;
}
