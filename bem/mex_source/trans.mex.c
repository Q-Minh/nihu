/*TRANS - perform translation step in the fmbem procedure
 *  N = TRANS(F, I, D, P, Perm, M)
 */

/* $Make: mex -O -output trans trans.mex.c fmbem.c vector.c$ */

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
    double *Fr, *Fi, *Mr, *Mi;
    int32_t *I, *D, *P, *Perm;
    int nclus, ns, nil, nm;
    double *Nr, *Ni;

    /* transfer input parameters */
    Fr = mxGetPr(prhs[0]);		/* Far field signature */
    Fi = mxGetPi(prhs[0]);
    nclus = mxGetN(prhs[0]);
    ns = mxGetM(prhs[0]);		/* number of directions (unit sphere) */
    I = (int32_t *)mxGetData(prhs[1]);		/* Interaction lists */
    nil = mxGetM(prhs[1]);
    D = (int32_t *)mxGetData(prhs[2]);		/* Distance indices lists */
    P = (int32_t *)mxGetData(prhs[3]);		/* Permutation indices lists */
    Perm = (int32_t *)mxGetData(prhs[4]);    /* Permutation vectors */
    Mr = mxGetPr(prhs[5]);		/* Translation operators */
    Mi = mxGetPi(prhs[5]);
    nm = mxGetN(prhs[5]);

    /* output vector */
    plhs[0] = mxCreateDoubleMatrix(ns, nclus, mxCOMPLEX);
    Nr = mxGetPr(plhs[0]);		/* Near field signature */
    Ni = mxGetPi(plhs[0]);

    /* call C subroutine */
    translate(nclus, ns, Fr, Fi, nil, I, D, P, Perm, nm, Mr, Mi, Nr, Ni);

    return;
}
