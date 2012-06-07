/*TRANS
 *
 *
 * Peter Fiala
 * 2009
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
    double *Fr, *Fi, *I, *D, *P, *Perm, *Mr, *Mi;
    int nclus, ns, nil, nm;
    double *Nr, *Ni;

    /* transfer input parameters */
    Fr = mxGetPr(prhs[0]);		/* Far field signature */
    Fi = mxGetPi(prhs[0]);
    nclus = mxGetN(prhs[0]);
    ns = mxGetM(prhs[0]);		/* number of directions (unit sphere) */
    I = mxGetPr(prhs[1]);		/* Interaction lists */
    nil = mxGetM(prhs[1]);
    D = mxGetPr(prhs[2]);		/* Distance indices lists */
    P = mxGetPr(prhs[3]);		/* Permutation indices lists */
    Perm = mxGetPr(prhs[4]);    /* Permutation vectors */
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

