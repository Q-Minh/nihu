/*bemHG_lin  Generate acoustic BEM system matrices with Burton Miller method
 *   [H, G] = bemHG_lin_bm(nodes, elements, g3, g4, dist, k, alpha)
 *   Computes the acoustic BEM system matrices H and G for the case of
 *   linear shape functions.
 * Parameters:
 *   nodes    : N x 3 matrix of xyz coordinates of model vertices
 *   elements : matrix of element node indices. Each row describes one
 *              element in the form [3 n1 n2 n3 0] for TRIA and
 *              [4 n1 n2 n3 n4] for QUAD elements.
 *   g3, g4   : Gaussian quadrature structure for TRIA and QUAD elements
 *              obtained from GAUSSQUAD2
 *   dist     : Distance limits governing integration density parameter
 *   k        : acoustic wave number
 *   alpha    : Burton-Miller coupling constant
 *
 */

/* $Make: mex -O -output bemHG_lin_bm bemHG_lin_bm.mex.c bemHG_li_bm.c integral_direct_bm.c mesh.c quadrature.c element.c vector.c green.c $ */

#include "mex.h"
#include "bemHG_li_bm.h"

/* ------------------------------------------------------------------------ */
/* The entry point of the mex function.                                     */
/* This function reads Matlab parameters, allocates space for the output    */
/* and calls C routines that perform the computations                       */
void mexFunction(int nlhs, 				/* number of lhs parameters */
                 mxArray *plhs[], 		/* pointer to lhs parameters */
                 int nrhs, 				/* number of rhs parameters */
                 const mxArray *prhs[]) /* pointer to rhs parameters */
{
    double *nodes, *elements, *points, *dist, k, alphar, alphai;
    int nnodes, nelements, npoints;
    gauss_t *g3, *g4;
    double *Ar, *Ai, *Br, *Bi;

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
        g3[j].Nxi = mxGetPr(mxGetField(prhs[2], j, "Nxi"));
        g3[j].Neta = mxGetPr(mxGetField(prhs[2], j, "Neta"));
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
    alphar = *mxGetPr(prhs[6]);
	alphai = *mxGetPi(prhs[6]);

	/* Allocate output parameters */
	plhs[0] = mxCreateDoubleMatrix(nnodes, nnodes, mxCOMPLEX);
	plhs[1] = mxCreateDoubleMatrix(nnodes, nnodes, mxCOMPLEX);
	Ar = mxGetPr(plhs[0]);
	Ai = mxGetPi(plhs[0]);
	Br = mxGetPr(plhs[1]);
	Bi = mxGetPi(plhs[1]);

	/* call C subroutine */
	/* NOTE: A: matrix H, B: matrix G */
	matrix_surf_lin_bm(nnodes, nodes, nelements, elements, g3, g4, dist, 
		k, alphar, alphai, Ar, Ai, Br, Bi);

    free(g3);
    free(g4);

    return;
}