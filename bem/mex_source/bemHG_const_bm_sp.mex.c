/*bemHG_const_bm_sp  Generate sparse acoustic BEM system matrices
 */

/* $Make: mex -O -output bemHG_const_bm_sp bemHG_const_bm_sp.mex.c bemHG_con_bm.c mesh.c integral_direct_bm.c quadrature.c element.c vector.c green.c$ */

#include "mex.h"
#include "bemHG_con_bm.h"
#include "types.h"

/* ------------------------------------------------------------------------ */
/* The entry point of the mex function. */
/* This function reads Matlab parameters, allocates space for the output */
/* and calls C routines that perform the computations */
void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
    /* input parameters */
    double *nodes, *elements, *points, *dist, k, alphar, alphai, *pairs;
    int nnodes, nelements, npoints, npairs;
    gauss_t *g3, *g4;
    /* output parameters */
    double *Ar, *Ai, *Br, *Bi;
    /* local variables */
    int j;

    /* Allocate memory for Gaussian integration structures */
    g3 = (gauss_t *)calloc(4, sizeof(gauss_t));
    g4 = (gauss_t *)calloc(4, sizeof(gauss_t));

    /* transfer input parameters */
    nodes = mxGetPr(prhs[0]);
    nnodes = mxGetM(prhs[0]);
    elements = mxGetPr(prhs[1]);
    nelements = mxGetM(prhs[1]);
    for (j = 0; j < 4; j++)
    {
        g3[j].num = mxGetM(mxGetField(prhs[2], j, "N"));
        g3[j].N = mxGetPr(mxGetField(prhs[2], j, "N"));
        g3[j].w = mxGetPr(mxGetField(prhs[2], j, "w"));
    }
    for (j = 0; j < 4; j++)
    {
        g4[j].num = mxGetM(mxGetField(prhs[3], j, "N"));
        g4[j].N = mxGetPr(mxGetField(prhs[3], j, "N"));
        g4[j].Nxi = mxGetPr(mxGetField(prhs[3], j, "Nxi"));
        g4[j].Neta = mxGetPr(mxGetField(prhs[3], j, "Neta"));
        g4[j].w = mxGetPr(mxGetField(prhs[3], j, "w"));
    }
    dist = mxGetPr(prhs[4]);
    k = mxGetScalar(prhs[5]);
	alphar = *(mxGetPr(prhs[6]));
	alphai = *(mxGetPi(prhs[6]));

	pairs = mxGetPr(prhs[7]);
	npairs = mxGetM(prhs[7]);

	/* Allocate output parameters */
	plhs[0] = mxCreateDoubleMatrix(npairs, 1, mxCOMPLEX);
	plhs[1] = mxCreateDoubleMatrix(npairs, 1, mxCOMPLEX);
	Ar = mxGetPr(plhs[0]);
	Ai = mxGetPi(plhs[0]);
	Br = mxGetPr(plhs[1]);
	Bi = mxGetPi(plhs[1]);

	/* call C subroutine */
	matrix_surf_const_bm_sparse(nnodes, nodes, nelements, elements,
							 npairs, pairs, g3, g4, dist, k, alphar, alphai, Ar, Ai, Br, Bi);

    free(g3);
    free(g4);

    return;
}
