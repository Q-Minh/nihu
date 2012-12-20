/*bemHG_const  Generate acoustic BEM system matrices
 *   [H, G] = bemHG_const(nodes, elements, g, dist, k)
 %   or
 *   [H, G] = bemHG_const(nodes, elements, g, dist, k, points)
 *   computes the acoustic BEM system matrices H and G for the case of
 *   constant shape functions.
 * Parameters:
 *   nodes    : N x 2 matrix of xy coordinates of model vertices
 *   elements : E x 2 matrix of element node indices
 *   g        : Gaussian quadrature structure for LINE elements
 *              obtained from GAUSSQUAD
 *   dist     : Distance limits governing integration density parameter
 *   k        : acoustic wave number
 *   points   : optional, if not given, then the acoustuc surface matrices
 *              are computed. If defined, then the field point matrices
 *              are returned. Points is a M x 2 matrix containing the xy
 *              coordinates of the field points
 *
 * Peter Fiala
 * 2012
 */

/* $Make: mex -O -output bemHG_const2D bemHG_const2D.mex.c bemHG_con.c mesh.c integral_direct.c element.c quadrature.c vector.c green.c$ */

#include "mex.h"
#include "bemHG_con.h"
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
    double *nodes, *elements, *points, *dist, k;
    int nnodes, nelements, npoints;
    gauss2D_t *g;
    /* output parameters */
    double *Ar, *Ai, *Br, *Bi;
    /* local variables */
    int j;

    /* Allocate memory for Gaussian integration structures */
    g = (gauss_t *)calloc(4, sizeof(gauss2D_t));

    /* transfer input parameters */
    nodes = mxGetPr(prhs[0]);
    nnodes = mxGetM(prhs[0]);
    elements = mxGetPr(prhs[1]);
    nelements = mxGetM(prhs[1]);
    for (j = 0; j < 4; j++)
    {
        g[j].num = mxGetM(mxGetField(prhs[2], j, "N"));
        g[j].N = mxGetPr(mxGetField(prhs[2], j, "N"));
        g[j].Nxi = mxGetPr(mxGetField(prhs[2], j, "Nxi"));
        g[j].w = mxGetPr(mxGetField(prhs[2], j, "w"));
    }
    dist = mxGetPr(prhs[3]);
    k = mxGetScalar(prhs[4]);

    if (nrhs == 5)
    {
        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(nelements, nelements, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(nelements, nelements, mxCOMPLEX);
        Ar = mxGetPr(plhs[0]);
        Ai = mxGetPi(plhs[0]);
        Br = mxGetPr(plhs[1]);
        Bi = mxGetPi(plhs[1]);

        /* call C subroutine */
        matrix_surf_const2D(nnodes, nodes, nelements, elements,
                          g, dist, k, Ar, Ai, Br, Bi);
    }
    else if (nrhs == 6)
    {
        points = mxGetPr(prhs[5]);
        npoints = mxGetM(prhs[5]);

        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(npoints, nelements, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(npoints, nelements, mxCOMPLEX);
        Ar = mxGetPr(plhs[0]);
        Ai = mxGetPi(plhs[0]);
        Br = mxGetPr(plhs[1]);
        Bi = mxGetPi(plhs[1]);

        /* call C subroutine */
        matrix_field_const2D(nnodes, nodes, nelements, elements,
                           npoints, points, g, dist, k, Ar, Ai, Br, Bi);
    }

    free(g);
    return;
}
