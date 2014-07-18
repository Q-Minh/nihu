/*leafHp
 *
 *
 * Peter Fiala
 * 2009
 */

/* $Make: mex -O -output leafHp leafHp.mex.c fmbem.c vector.c $ */

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
    double *r, *n, *pr, *pi, *R, *s, k;
    int32_t *father;
    int nnod, nclus, ns;
    double *Fr, *Fi;

    /* transfer input parameters */
    r = mxGetPr(prhs[0]);		/* source nodes */
    nnod = mxGetN(prhs[0]);
    n = mxGetPr(prhs[1]);		/* source normals */
    pr = mxGetPr(prhs[2]);		/* excitation pressure */
    pi = mxGetPi(prhs[2]);
    father = (int32_t *)mxGetData(prhs[3]);	/* father cluster indices */
    nclus = (int)mxGetScalar(prhs[4]);		/* cluster centres */
    s = mxGetPr(prhs[5]);		/* directions (unit sphere) */
    ns = mxGetN(prhs[5]);
    k = mxGetScalar(prhs[6]);	/* wave number */

    /* Allocate output parameters */
    plhs[0] = mxCreateDoubleMatrix(ns, nclus, mxCOMPLEX);
    Fr = mxGetPr(plhs[0]);		/* Far field signature */
    Fi = mxGetPi(plhs[0]);

    /* call C subroutine */
    leafHp(nnod, r, n, pr, pi, father, nclus, ns, s, k, Fr, Fi);

    return;
}
